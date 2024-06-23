# Device objects
# Sam Chinnery
# 2022-01-27

# Contains material properties and carrier distribution functions
mutable struct Model

	e0::Float64		# Vacuum permittivity (constant)
	q::Float64		# Electron charge (constant)
	vt::Float64		# Thermal voltage at 300 K (constant)
	er::Dict{Int64,Float64} # Relative permittivity by cell region
	ni::Dict{Int64,Float64}	# Intrinsic carrier concentration by cell region
	tn::Dict{Int64,Float64}	# Electron lifetime by cell region
	tp::Dict{Int64,Float64}	# Hole lifetime by cell region
	ni_boundary::Dict{Int64,Float64}
				# Intrinsic concentration by boundary region
	ew::Dict{Int64,Float64}	# Metal-insulator workfunction by cell region

	doping			# Doping profile model (function)
	recomb			# Recombination model (function)
	generation		# Generation model (function)
	mobility_n		# Electron mobility model (function)
	mobility_p		# Hole mobility model (function)
	n_dist			# Electron statistical model (function)
	p_dist			# Hole statistical model (function)
	doping_to_np		# Doping to charge-neutral n, p (function)
	np_to_v			# Charge-neutral n, p to potential (function)

	fldmob::Bool		# Enable/disable field-dependent mobility
	impact::Bool		# Enable/disable impact ionization

end

# Default constructor with empty properties
Model() = Model(
	8.854_187_812_8e-18,
	1.602_176_634e-19,
	25.851_999_786_435_5e-3,
	props_default,
	props_default,
	props_default,
	props_default,
	props_default,
	props_default,
	doping_default,
	u_srh,
	g_impact_selbherr,
	un_vsat,
	up_vsat,
	n_boltzmann,
	p_boltzmann,
	np_bc,
	v_bc,
	false,
	false
)

# Contains device model, geometry and maps of region numbers to region types
struct Device
	model::Model				# Device model and parameters
	grid::ExtendableGrids.ExtendableGrid	# FVM discretization grid
	b_types::Dict{Int64,String}		# Boundary region types
	r_types::Dict{Int64,String}		# Cell region types
end

# Default initializers for Model
props_default = Dict(0=>0.0)
doping_default(x,y,reg) = 0
u_srh(n,p,ni,tn,tp) = (n*p-ni^2)/((n+ni)*tp+(p+ni)*tn)
n_boltzmann(v,qfp_n,ni,vt) = ni*exp((v-qfp_n)/vt)
p_boltzmann(v,qfp_p,ni,vt) = ni*exp((qfp_p-v)/vt)
v_bc(n0,p0,ni,vt) = n0>p0 ? vt*log(n0/ni) : -vt*log(p0/ni)

# Compute equilibrium concentrations using the charge neutrality assumption for
# Ohmic contacts. Larger quadratic root is used to avoid instability
function np_bc(gamma,ni)
	if gamma>0.0
		n0 = 1/2*(gamma+sqrt(gamma^2+4*ni^2))
		p0 = ni^2/n0
	elseif gamma<0.0
		p0 = 1/2*(-gamma+sqrt(gamma^2+4*ni^2))
		n0 = ni^2/p0
	else
		n0 = ni
		p0 = ni
	end
	return (n0,p0)
end

# Velocity-saturation mobility models of Caughey and Thomas (1967), "Carrier
# mobilities in silicon empirically related to doping and field"
function un_vsat(n,p,e_norm,reg,fldmob)
	u0 = 1.375e11
	if !fldmob
		return (u0)
	else
		e_crit = 1.95e0
		beta = 2.0
		u = u0/(1+(e_norm/e_crit)^beta)^(1/beta)
		return (u)
	end
end
function up_vsat(n,p,e_norm,reg,fldmob)
	u0 = 4.87e10
	if !fldmob
		return (u0)
	else
		e_crit = 8e-1
		beta = 1.0
		u = u0/(1+(e_norm/e_crit)^beta)^(1/beta)
		return (u)
	end
end

# Impact ionization model of Thornber (1981), "Applications of scaling to
# problems in high-field electronic transport"
function g_impact_thornber(e_norm,jn_norm,jp_norm,device)

	# Extract model from device
	m = device.model

	# Using the parameter set of Van Overstraeten and DeMan
	ei_n = 3.6
	ei_p = 6.2
	fi_n = 1.404e2
	fi_p = 5.933e2
	fr_n = 2.229e1
	fr_p = 5.910e0
	fkt_n = 9.747e-1
	fkt_p = 2.392e0

	# Ionization coefficients
	alpha_n = e_norm/ei_n * exp(-fi_n/(e_norm*(1+e_norm/fr_n)+fkt_n))
	alpha_p = e_norm/ei_p * exp(-fi_p/(e_norm*(1+e_norm/fr_p)+fkt_p))

	# Ionization rate
	g = 1/m.q * (alpha_n*jn_norm+alpha_p*jp_norm)
	return (g)

# function g_impact_thornber
end

# Impact ionization model of Selbherr (1986), "On modeling MOS-devices"
function g_impact_selbherr(
	e_norm,
	jn_norm,
	jp_norm,
	device
)

	# Extract model from device
	m = device.model

	# From reference [26] in Selbherr
	an_inf = 1.0e2
	ap_inf = 2.0e2
	ecrit_n = 1.66e2
	ecrit_p = 1.98e2
	gamma_n = 1.0
	gamma_p = 1.0

	# Ionization coefficients
	alpha_n = an_inf*exp(-ecrit_n/e_norm*gamma_n)
	alpha_p = ap_inf*exp(-ecrit_p/e_norm*gamma_p)

	# Ionization rate
	g = 1/m.q * (alpha_n*jn_norm+alpha_p*jp_norm)
	return (g)

# function g_impact_selbherr
end

