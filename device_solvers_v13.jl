# Equilibrium and non-equilibrium solvers for Semiconductors.jl
# Sam Chinnery
# 2022-01-28

# Charge neutrality initial conditions
function ic_equilib(d,sys::VoronoiFVM.AbstractSystem{Tv,Ti}) where {Tv,Ti}

	# Extract model information and create initial condition vector
	m = d.model
	ic = unknowns(sys)
	ic .= 0.0

	# Get geometry data
	dims = dim_space(d.grid)
	node = VoronoiFVM.Node{Tv,Ti}(sys)
	geom = sys.grid[CellGeometries][1]
	nodes_per_cell = num_nodes(geom)
	n_cells = num_cells(sys.grid)

	# Loop through each discretization cell
	for icell in 1: n_cells

		# Get parameters that are constant over each cell
		reg = sys.grid[CellRegions][icell]
		if d.r_types[reg]=="semiconductor"
			ni = m.ni[reg]
		end

	# Loop through each node in the current cell
	for inode in 1: nodes_per_cell

		# Populate node object with adjacent node information and index
		# in solution vector
		VoronoiFVM._fill!(node,inode,icell)
		inode_k = node.index

		# Compute equilibrium potential assuming charge neutrality
		if d.r_types[reg]=="semiconductor"
			doping = m.doping(node[1],node[2],reg,boundary=false)
			n0,p0 = m.doping_to_np(doping,ni)
			v = m.np_to_v(n0,p0,ni,m.vt)
		else
			v = 0.0
		end

		# Update the initial condition, or average it with whatever is
		# already there
		if ic[1,node.index] == 0.0
			ic[1,node.index] = v
		else
			#ic[1,node.index] = v
			ic[1,node.index] = (ic[1,node.index]+v)/2.0
		end

	# for inode in 1: nodes_per_cell
	end

	# for icell in 1: n_cells
	end

	return (ic)

# function ic_equilib
end

# Equilibrium solver
function equilib(device;Plotter=nothing,verbose=false,unknown_storage=:sparse,
	damp_initial=0.1)

	# Extract model and grid from Device object
	m = device.model
	grid = device.grid

	# Define flux, source and reaction terms for Poisson-Boltzmann equation
	function flux!(f,u,edge)
		eps = m.e0*m.er[edge.region]
		f[1] = eps*(u[1,1]-u[1,2])
	end
	function source!(f,node)
		if device.r_types[node.region]=="semiconductor"
			f[1] = m.q * m.doping(node[1],node[2],node.region)
		end
	end
	function reaction!(f,u,node)
		if device.r_types[node.region]=="semiconductor"
			ni = m.ni[node.region]
			f[1] = 2*m.q*ni * sinh(u[1]/m.vt)
		end
	end

	# Dirichlet condition for potential at each boundary node
	function bcond!(f,u,bnode)

		# Region of boundary node
		reg = bnode.region

		# Neumann BCs are assumed at non-contact boundaries
		if device.b_types[reg]=="contact"

			# Compute equilibrium carrier densities for each
			# boundary node
			ni = m.ni_boundary[reg]
			n0,p0 = m.doping_to_np(m.doping(
				bnode[1],bnode[2],reg,boundary=true),ni)
			v = m.np_to_v(n0,p0,ni,m.vt)

			# Apply Dirichlet BCs for potential and concentrations
			boundary_dirichlet!(f,u,bnode,1,reg,v)

		# At gate contacts, potential is constrained by the
		# metal-insulator workfunction
		elseif device.b_types[reg]=="gate"

			# Look up workfunction and set Dirichlet condition for
			# potential only
			boundary_dirichlet!(f,u,bnode,1,reg,m.ew[reg])

		end

	end

	# Create a finite volume system
	sys = VoronoiFVM.System(grid,flux=flux!,source=source!,
		reaction=reaction!,bcondition=bcond!,species=[1],
		unknown_storage=unknown_storage)

	# Use electroneutrality assumption everywhere to compute IC
	inival = ic_equilib(device,sys)
	solution = unknowns(sys)

	# Solve system
	control = VoronoiFVM.SolverControl()
	control.verbose = verbose
	control.damp_initial = damp_initial
	solve!(solution,inival,sys,control=control)

	# Plot solution
	scalarplot(grid,solution[1,: ],title="Equilibrium potential", 
		Plotter=Plotter,colormap=:coolwarm,resolution=(1000,500))

	return (solution)

end

# Physics for non-equilibrium solver, using quasi-Fermi potentials
# TODO: This is currently broken
function physics_noneq_qfp(device)

	# Extract model and grid from Device object
	m = device.model
	grid = device.grid

	# Flux between control volumes, using quasi-Fermi potentials
	function flux!(f,u,edge)

		# Poisson flux
		eps = m.e0*m.er[edge.region]
		f[1] = eps*(u[1,1]-u[1,2])

		if device.r_types[edge.region]=="semiconductor"

			# Bernoulli evaluations for S-G
			bp,bm = fbernoulli_pm((u[1,2]-u[1,1])/m.vt)

			# Compute carrier concentrations from QFPs
			ni = m.ni[edge.region]
			fn_k = m.n_dist(u[1,1],u[2,1],ni,m.vt)
			fn_l = m.n_dist(u[1,2],u[2,2],ni,m.vt)
			fp_k = m.p_dist(u[1,1],u[3,1],ni,m.vt)
			fp_l = m.p_dist(u[1,2],u[3,2],ni,m.vt)

			# QFP flux equations
			un = m.un[edge.region]
			up = m.up[edge.region]
			f[2] = un*m.vt * ( fn_k*bm - fn_l*bp )
			f[3] = un*m.vt * ( fp_k*bp - fp_l*bm )

		end

	end

	# Source term, only nonzero in Poisson equation
	function source!(f,node)
		if device.r_types[node.region]=="semiconductor"
			f[1] = m.q*m.doping(node[1],node[2],node.region)
		end
	end

	# Reaction term, using quasi-Fermi potentials
	function reaction!(f,u,node)
		if device.r_types[node.region]=="semiconductor"
			ni = m.ni[node.region]
			n_val = m.n_dist(u[1],u[2],ni,m.vt)
			p_val = m.n_dist(u[1],u[3],ni,m.vt)
			tn = m.tn[node.region]
			tp = m.tp[node.region]
			recomb = m.recomb(n_val,p_val,ni,tn,tp)
			f[1] = m.q*(n_val-p_val)
			f[2] = recomb
			f[3] = recomb
		end
	end

	# Boundary conditions, using quasi-Fermi potentials
	function bcond!(f,u,bnode,bias)

		# Region of boundary node
		reg = bnode.region

		# Neumann BCs are assumed at non-contact boundaries
		if device.b_types[reg] == "contact"

			# Compute equilibrium carrier densities for each
			# boundary node
			ni = m.ni_boundary[reg]
			n0,p0 = m.doping_to_np(m.doping(
				bnode[1],bnode[2],reg,boundary=true),ni)
			v = m.np_to_v(n0,p0,ni,m.vt)

			# Apply Dirichlet BCs for potential and concentrations
			boundary_dirichlet!(f,u,bnode,1,reg,v+bias[reg])
			boundary_dirichlet!(f,u,bnode,2,reg,bias[reg])
			boundary_dirichlet!(f,u,bnode,3,reg,bias[reg])
		end

	end

	# Return physics struct
	return (flux!,source!,reaction!,bcond!)

end

# Physics for non-equilibrium solver, using Boltzmann statistics for carrier
# concentrations
function physics_noneq_boltz(device)

	# Extract model and grid from Device object
	m = device.model
	grid = device.grid

	# Flux between control volumes, using Boltzmann statistics for
	# carrier concentrations
	function flux!(f,u,edge)

		# This will only work if field-dependent mobility is not enabled
		@assert !m.fldmob

		# Poisson flux
		eps = m.e0*m.er[edge.region]
		f[1] = eps*(u[1,1]-u[1,2])

		# Carrier concentration flux equations
		if device.r_types[edge.region]=="semiconductor"
			bp,bm = fbernoulli_pm((u[1,2]-u[1,1])/m.vt)
			un = m.mobility_n(0,0,0,edge.region,m.fldmob)
			up = m.mobility_p(0,0,0,edge.region,m.fldmob)

			#sigma(x) = 1.0/(1.0+exp(0.5*x))
			#dv = (u[1,2]-u[1,1])/m.vt
			#bm = 1.0+dv*(1.0-sigma(dv))
			#bp = 1.0-dv*sigma(dv)

			f[2] = un*m.vt * ( u[2,1]*bm - u[2,2]*bp )
			f[3] = up*m.vt * ( u[3,1]*bp - u[3,2]*bm )
		end

	end

	# Source term, only nonzero in Poisson equation
	function source!(f,node)
		if device.r_types[node.region]=="semiconductor"
			f[1] = m.q*m.doping(node[1],node[2],node.region)
		end
	end

	# Reaction term, using Boltzmann statistics for carrier concentrations
	function reaction!(f,u,node)
		if device.r_types[node.region]=="semiconductor"
			ni = m.ni[node.region]
			tn = m.tn[node.region]
			tp = m.tp[node.region]
			recomb = m.recomb(u[2],u[3],ni,tn,tp)
			f[1] = m.q*(u[2]-u[3])
			f[2] = recomb
			f[3] = recomb
		end
	end

	# Boundary conditions, using carrier concentrations
	function bcond!(f,u,bnode)

		# Region of boundary node
		reg = bnode.region

		# Neumann BCs are assumed at non-contact boundaries
		if device.b_types[reg]=="contact"

			# Bias voltage of current contact
			bias = parameters(u)[reg]

			# Compute equilibrium carrier densities for each
			# boundary node
			ni = m.ni_boundary[reg]
			n0,p0 = m.doping_to_np(m.doping(
				bnode[1],bnode[2],reg,boundary=true),ni)
			v = m.np_to_v(n0,p0,ni,m.vt)

			# Apply Dirichlet BCs for potential and concentrations
			boundary_dirichlet!(f,u,bnode,1,reg,v+bias)
			boundary_dirichlet!(f,u,bnode,2,reg,n0)
			boundary_dirichlet!(f,u,bnode,3,reg,p0)

		# At gate contacts, potential is constrained by the
		# metal-insulator workfunction
		elseif device.b_types[reg]=="gate"

			# Bias voltage of current contact
			bias = parameters(u)[reg]

			# Look up workfunction and set Dirichlet condition for
			# potential only
			boundary_dirichlet!(f,u,bnode,1,reg,m.ew[reg]+bias)

		end

	end

	# Return physics functions
	return (flux!,source!,reaction!,bcond!)

end

# Compute norm of the electric field, assumed constant over each grid cell
function cell_field_norms(u,sys,device,idx,icell,edge,edges_per_cell)

	# Extract model from device
	m = device.model

	# Matrices for computing field magnitudes
	dims = dim_space(sys.grid)
	tangents = zeros(Float64,edges_per_cell,dims)
	e_proj = zeros(typeof(u[1]),edges_per_cell,1)
	jn_proj = zeros(typeof(u[1]),edges_per_cell,1)
	jp_proj = zeros(typeof(u[1]),edges_per_cell,1)

	# Electric field loop
	for iedge in 1: edges_per_cell

		# Get adjacent nodes for each edge
		VoronoiFVM._fill!(edge,iedge,icell)
		inode_k = edge.node[1]
		inode_l = edge.node[2]

		# Compute length of edge and tangent unit vector
		xk = sys.grid[Coordinates][: ,inode_k]
		xl = sys.grid[Coordinates][: ,inode_l]
		hkl = norm(xl-xk)
		tangent = (xl-xk)./hkl

		# Indices of V in the solution vector
		idx_vk = idx[1,inode_k]
		idx_vl = idx[1,inode_l]

		# Populate tangent and projection matrices
		tangents[iedge,: ] = tangent
		e_proj[iedge] = (u[idx_vl]-u[idx_vk])/hkl

	end

	# Compute E field and its magnitude using least-squares solution (since
	# the tangent matrix is tall rectangular)
	e_field = tangents\e_proj
	e_norm = norm(e_field)

	# Only do currents if impact ionization enabled
	jn_norm = 0.0
	jp_norm = 0.0
	if m.impact

	# Electron and hole current loop
	for iedge in 1: edges_per_cell

		# Get adjacent nodes for each edge
		VoronoiFVM._fill!(edge,iedge,icell)
		inode_k = edge.node[1]
		inode_l = edge.node[2]

		# Compute length of edge and tangent unit vector
		xk = sys.grid[Coordinates][: ,inode_k]
		xl = sys.grid[Coordinates][: ,inode_l]
		hkl = norm(xl-xk)

		# Indices of V, n and p in the solution vector
		idx_vk = idx[1,inode_k]
		idx_vl = idx[1,inode_l]
		idx_nk = idx[2,inode_k]
		idx_nl = idx[2,inode_l]
		idx_pk = idx[3,inode_k]
		idx_pl = idx[3,inode_l]

		# Compute intermediate values for current discretization
		dv = (u[idx_vl]-u[idx_vk])/m.vt
		bp,bm = fbernoulli_pm(dv)
		un = m.mobility_n(0,0,e_norm,edge.region,m.fldmob)
		up = m.mobility_p(0,0,e_norm,edge.region,m.fldmob)

		# Populate projection matrices
		jn_proj[iedge] = m.q*un*m.vt/hkl*(u[idx_nk]*bm-u[idx_nl]*bp)
		jp_proj[iedge] = m.q*up*m.vt/hkl*(u[idx_pk]*bp-u[idx_pl]*bm)

	end

	# Compute J fields and their magnitudes
	jn_field = tangents\jn_proj
	jp_field = tangents\jp_proj
	jn_norm = norm(jn_field)
	jp_norm = norm(jp_field)

	# if impact
	end

	# Compute and show residual (debugging only)
	#resid = norm(tangents*e_field-e_proj)/e_norm
	#print(value(resid)," ")

	return (e_norm,jn_norm,jp_norm)

end

# C(x) = (1-B(x)) / x = 1/x - 1/(e^x-1)
# Also referred to as Q(x) in Laux and Byrnes (1985), "Semiconductor simulation
# using generalized mobility models"
function c_horner(x)
	y = x/47_900_160
	y = -x*y
	y = x*(1/1_209_600+y)
	y = -x*y
	y = x*(1/30_240+y)
	y = -x*y
	y = x*(1/720+y)
	y = -x*y
	y = x*(1/12+y)
	y = 1/2-y
end

# Positive and negative evaluations of C(x)
function c_pm(x)
	large_threshold = 40.0
	small_threshold = 0.25
	if x<-large_threshold
		return (1.0+1.0/x,-1.0/x)
	elseif x>large_threshold
		return (1.0/x,1.0-1.0/x)
	else
		expx = exp(x)
		denom = 1.0-expx
		if abs(x)>small_threshold
			bp = 1.0/denom+1.0/x
			return (bp,1.0-bp)
		else
			bp = c_horner(x)
			return (bp,1.0-bp)
		end
	end
end

# Bernoulli function needs to handle symbolic arguments to allow automatic
# sparsity detection
fbernoulli_pm(x) = VoronoiFVM.fbernoulli_pm(x)
fbernoulli_pm(x::Symbolics.Num) = (x/(exp(x)-1),-x/(exp(-x)-1))
c_pm(x::Symbolics.Num) = (1/x-1/(exp(x)-1),-1/x-1/(exp(-x)-1))

# Generic operator for discretizing flux with field-dependent mobility
# Will compute boundary integral if a test function is given in tf
function _generic_noneq!(f,u::AbstractArray{Tu},
	sys::VoronoiFVM.AbstractSystem{Tv,Ti},device;
	tf=nothing) where {Tu,Tv,Ti}

	# Extract semiconductor model data
	m = device.model

	# Get unknown indices
	idx = unknown_indices(unknowns(sys))

	# If test function not given, initialize residual vector. Otherwise, set
	# up for integration
	if tf==nothing
		f .= 0.0
	else

		# Holds integral when test function is given. Poisson part of
		# boundary integral is set to 0 since we don't need it
		integral = zeros(Tu,3)

		# We need to reshape U for the integration so that NaNs are
		# dropped in the insulator region and idx can be used for
		# indexing
		u = VoronoiFVM.values(u)

	end

	# Get grid geometry, total number of grid cells
	# (segments/triangles/tets) and number of edges per grid cell
	geom = sys.grid[CellGeometries][1]
	n_cells = num_cells(sys.grid)
	nodes_per_cell = num_nodes(geom)
	edges_per_cell = num_edges(geom)

	# Create node and edge objects to hold adjacent node information
	node = VoronoiFVM.Node{Tv,Ti}(sys)
	edge = VoronoiFVM.Edge{Tv,Ti}(sys)

	# Loop through each grid cell
	for icell in 1: n_cells

	# Material parameters that are constant in each grid cell
	region = sys.grid[CellRegions][icell]
	eps = m.e0*m.er[region]

	# Only compute E field if field-dependent mobility is enabled
	if device.r_types[region]=="semiconductor"
		e_norm,jn_norm,jp_norm = cell_field_norms(
			u,sys,device,idx,icell,edge,edges_per_cell)
	else
		e_norm = 0.0
		jn_norm = 0.0
		jp_norm = 0.0
	end

	# Only need reaction part if impact ionization is enabled
	if m.impact

	# Reaction discretization loop
	for inode in 1: nodes_per_cell

		# Populate node object with the current node
		VoronoiFVM._fill!(node,inode,icell)
		inode_k = node.index
		node_factor = sys.cellnodefactors[inode,icell]

		# Index of V, n, p in the solution vector
		idx_vk = idx[1,inode_k]
		idx_nk = idx[2,inode_k]
		idx_pk = idx[3,inode_k]

		# Electron and hole reaction terms (only in semiconductor
		# regions)
		if device.r_types[region]=="semiconductor"

			# Poisson reaction term
			v_res = m.q*(u[idx_nk]-u[idx_pk])

			# Compute recombination
			ni = m.ni[node.region]
			tn = m.tn[node.region]
			tp = m.tp[node.region]
			recomb = m.recomb(u[idx_nk],u[idx_pk],ni,tn,tp)
			generation = m.generation(e_norm,jn_norm,jp_norm,device)

			# Add reaction residual if test function not given; add
			# curent part of integral if it is
			if tf==nothing
				carrier_res = node_factor*(recomb-generation)
				f[idx_vk] += node_factor*v_res
				f[idx_nk] += carrier_res
				f[idx_pk] += carrier_res
			else
				carrier_int = node_factor*(recomb-generation)*
					tf[inode_k]
				integral[2] += carrier_int
				integral[3] += carrier_int
			end

		# if device.r_types[region]=="semiconductor"
		end

	# for inode in 1: nodes_per_cell
	end

	# if m.impact
	end

	# Only do flux discretization if field-dependent mobility is enabled
	if m.fldmob

	# Flux discretization loop
	for iedge in 1: edges_per_cell

		# Populate edge object with adjacent nodes
		VoronoiFVM._fill!(edge,iedge,icell)

		# Get indices of adjacent nodes, and get the factor sigma/h for
		# flux discretization
		inode_k = edge.node[1]
		inode_l = edge.node[2]
		edge_factor = sys.celledgefactors[iedge,icell]

		# Indices of V in the solution vector
		idx_vk = idx[1,inode_k]
		idx_vl = idx[1,inode_l]

		# Poisson flux discretization
		if tf==nothing
			v_flux = eps*(u[idx_vk]-u[idx_vl])
			f[idx_vk] += edge_factor*v_flux
			f[idx_vl] -= edge_factor*v_flux
		end

		# Electron and hole flux discretization (only in semiconductor
		# regions)
		if device.r_types[region]=="semiconductor"

			# Indices of n, p in the solution vector
			idx_nk = idx[2,inode_k]
			idx_nl = idx[2,inode_l]
			idx_pk = idx[3,inode_k]
			idx_pl = idx[3,inode_l]

			# TODO: Use average carrier concentrations along
			# edge with C function
			un = m.mobility_n(0,0,e_norm,region,m.fldmob)
			up = m.mobility_p(0,0,e_norm,region,m.fldmob)

			# Continuity equations
			bp,bm = fbernoulli_pm((u[idx_vl]-u[idx_vk])/
				m.vt)
			n_flux = un*m.vt*(u[idx_nk]*bm-u[idx_nl]*bp)
			p_flux = up*m.vt*(u[idx_pk]*bp-u[idx_pl]*bm)

			# Add flux residual if test function not given; add
			# curent part of integral if it is
			if tf==nothing
				f[idx_nk] += edge_factor*n_flux
				f[idx_nl] -= edge_factor*n_flux
				f[idx_pk] += edge_factor*p_flux
				f[idx_pl] -= edge_factor*p_flux
			else
				tf_flux = tf[inode_k]-tf[inode_l]
				integral[2] += edge_factor*n_flux*tf_flux
				integral[3] += edge_factor*p_flux*tf_flux
			end

		# if device.r_types[region]=="semiconductor"
		end

	# for iedge in 1: edges_per_cell
	end

	# if m.fldmob
	end

	# for icell in 1: n_cells
	end

	# Return integral if test function is given
	if tf==nothing
		return (nothing)
	else
		return (integral)
	end

end

# Generic operator sparsity pattern
function _generic_noneq_sparsity(sys::VoronoiFVM.AbstractSystem{Tv,Ti},
	device) where {Tv,Ti}

	# Extract model from device
	m = device.model

	# Create empty sparse matrix to hold sparsity pattern
	sparsity = spzeros(num_dof(sys),num_dof(sys))

	# Get unknown indices
	idx = unknown_indices(unknowns(sys))

	# Get grid geometry, total number of grid cells
	# (segments/triangles/tets) and number of edges per grid cell
	geom = sys.grid[CellGeometries][1]
	n_cells = num_cells(sys.grid)
	nodes_per_cell = num_nodes(geom)
	edges_per_cell = num_edges(geom)

	# Create node and edge objects to hold adjacent node information
	node = VoronoiFVM.Node{Tv,Ti}(sys)
	edge = VoronoiFVM.Edge{Tv,Ti}(sys)

	# Loop through each cell
	for icell in 1: n_cells
	region = sys.grid[CellRegions][icell]

	# Only need reaction part if impact ionization is enabled
	if m.impact

	# Loop for reaction operator sparsity
	for inode in 1: nodes_per_cell

		# Populate node object with the current node
		VoronoiFVM._fill!(node,inode,icell)
		inode_k = node.index
		node_factor = sys.cellnodefactors[inode,icell]

		# Index of V, n, p in the solution vector
		idx_vk = idx[1,inode_k]
		idx_nk = idx[2,inode_k]
		idx_pk = idx[3,inode_k]

		# Reaction term sparsity, only in semiconductor regions
		if device.r_types[region]=="semiconductor"
			sparsity[idx_vk,idx_nk] = 1
			sparsity[idx_vk,idx_pk] = 1
			sparsity[idx_nk,idx_nk] = 1
			sparsity[idx_nk,idx_pk] = 1
			sparsity[idx_pk,idx_nk] = 1
			sparsity[idx_pk,idx_pk] = 1
		end

	# for inode in 1: nodes_per_cell
	end

	# if m.impact
	end

	# Only do flux discretization if field-dependent mobility is enabled
	if m.fldmob

	# Loop for flux operator sparsity
	for iedge in 1: edges_per_cell

		# Populate edge object with adjacent nodes
		VoronoiFVM._fill!(edge,iedge,icell)

		# Get indices of adjacent nodes
		inode_k = edge.node[1]
		inode_l = edge.node[2]

		# Indices of V in the solution vector
		idx_vk = idx[1,inode_k]
		idx_vl = idx[1,inode_l]

		# Poisson flux sparsity
		sparsity[idx_vk,idx_vk] = 1
		sparsity[idx_vk,idx_vl] = 1
		sparsity[idx_vl,idx_vl] = 1
		sparsity[idx_vl,idx_vk] = 1

		# Only define electron and hole sparsity in semiconductor
		# regions
		has_nk = sys.node_dof[2,inode_k]==2
		has_nl = sys.node_dof[2,inode_l]==2
		if device.r_types[region]=="semiconductor"

			# Indices of n, p in the solution vector
			idx_nk = idx[2,inode_k]
			idx_nl = idx[2,inode_l]
			idx_pk = idx[3,inode_k]
			idx_pl = idx[3,inode_l]

			# Electron flux sparsity
			sparsity[idx_nk,idx_nk] = 1
			sparsity[idx_nk,idx_nl] = 1
			sparsity[idx_nk,idx_vk] = 1
			sparsity[idx_nk,idx_vl] = 1
			sparsity[idx_nl,idx_nl] = 1
			sparsity[idx_nl,idx_nk] = 1
			sparsity[idx_nl,idx_vk] = 1
			sparsity[idx_nl,idx_vl] = 1

			# Hole flux sparsity
			sparsity[idx_pk,idx_pk] = 1
			sparsity[idx_pk,idx_pl] = 1
			sparsity[idx_pk,idx_vk] = 1
			sparsity[idx_pk,idx_vl] = 1
			sparsity[idx_pl,idx_pl] = 1
			sparsity[idx_pl,idx_pk] = 1
			sparsity[idx_pl,idx_vk] = 1
			sparsity[idx_pl,idx_vl] = 1

		# if device.r_types[region]=="semiconductor"
		end

	# for iedge in 1: edges_per_cell
	end

	# if m.fldmob
	end

	# for icell in 1: n_cells
	end

	# Plot the sparsity pattern for debugging purposes
	#matshow(sparsity[1: 9,1: 9],cmap=:binary)
	#xticks(0: 8,1: 9)
	#yticks(0: 8,1: 9)

	return (sparsity)

end

# Copied from VoronoiFVM.jl to allow plot of sparsity pattern for debugging
function _auto_sparsity(sys,g)
	generic_operator(f,u) = g(f,u,sys)
	input = rand(VoronoiFVM.num_dof(sys))
	output = similar(input)
	tdetect = @elapsed begin
		sparsity_pattern = Symbolics.jacobian_sparsity(
			generic_operator,output,input)
		sparsity = Float64.(sparse(sparsity_pattern))
	end
	println("sparsity detection for generic operator: $(tdetect) s")
	if nnz(sparsity)==0
		error("Sparsity detection failed: no pattern found")
	end
	#n = 48
	#matshow(sparsity[1: n,1: n],cmap=:binary)
	#xticks(0: n-1,1: n)
	#yticks(0: n-1,1: n)
	return (sparsity)
end

# Concentration-weighted test function discretization of Nanz et al. 1992
# ("Calculation of Contact Currents in Device Simulation")
function tf_sys(system::VoronoiFVM.AbstractSystem{Tv},sol,device,bc0,bc1;
	jp=false) where Tv

	# Extract variables from model and solution
	m = device.model
	v = sol[1,: ]
	n = sol[2,: ]
	p = sol[3,: ]

	# Flux discretization using S-G interpolation for concentrations along
	# edge
	function flux!(f,u,edge)

		# Node indices
		inode1 = edge.node[1]
		inode2 = edge.node[2]

		# C function evaluation
		dv = (v[inode2]-v[inode1])/m.vt
		cp,cm = c_pm(dv)

		# Flux discretization
		if !jp
			nk = n[inode1]
			nl = n[inode2]
			f[1] =  (cm*nk+cp*nl)*(u[2]-u[1])/m.ni[edge.region]
		else
			pk = p[inode1]
			pl = p[inode2]
			f[1] =  (cp*pk+cm*pl)*(u[2]-u[1])/m.ni[edge.region]
		end

	end

	# Apply 1.1 BCs at active contact and -0.1 at other contacts
	function bcond!(f,u,bnode)
		for ireg in bc1
			boundary_dirichlet!(f,u,bnode,1,ireg,1.1)
		end
		for ireg in bc0
			boundary_dirichlet!(f,u,bnode,1,ireg,-0.1)
		end
	end

	# Return test function system
	tfsystem = VoronoiFVM.System(system.grid,flux=flux!,bcondition=bcond!,
		species=[1],unknown_storage=:dense)
	return (VoronoiFVM.TestFunctionFactory{Tv}(system,tfsystem))

end

# Concentration-weighted test function solution
function tf_solve(factory::VoronoiFVM.TestFunctionFactory{Tv},ic) where Tv

	# Create initial condition and solution vectors
	inival = unknowns(factory.tfsystem)
	inival .= ic
	solution = unknowns(factory.tfsystem)

	# Solve test function system and apply clipping and sin^2 smoothing
	solve!(solution,inival,factory.tfsystem)
	sol_clipped = min.(max.(solution,0),1)
	#sol_smooth = sin.(pi/2*sol_clipped).^2
	sol_smooth = -2*sol_clipped.^3 .+ 3*sol_clipped.^2
	return vec(sol_smooth)

end

# Copied from VoronoiFVM.jl to allow use of generic operator
function integrate_generic(system::VoronoiFVM.AbstractSystem{Tv,Ti},
	tf::Vector{Tv},U::AbstractArray{Tu,2},device,idx) where {Tu,Tv,Ti}

	grid = system.grid
	nspecies = VoronoiFVM.num_species(system)
	res = zeros(Tu,nspecies)
	src = zeros(Tu,nspecies)

 	physics = system.physics
	node = VoronoiFVM.Node{Tv,Ti}(system)
	bnode = VoronoiFVM.BNode{Tv,Ti}(system)
	edge = VoronoiFVM.Edge{Tv,Ti}(system)
	bedge = VoronoiFVM.BEdge{Tv,Ti}(system)
	VoronoiFVM.@create_physics_wrappers(physics,node,bnode,edge,bedge)

	UK = Array{Tu,1}(undef,nspecies)
	UKL = Array{Tu,1}(undef,2*nspecies)
	geom = grid[CellGeometries][1]

	# Compute flux using generic operator
	integral = _generic_noneq!(nothing,U,system,device,tf=tf)

	# TODO: We loop through all species here, but it is not necessary to
	# integrate the Poisson equation to get the current
	for icell in 1: num_cells(grid)

	# Add contribution of source and reaction term
	for inode in 1: num_nodes(geom)
		VoronoiFVM._fill!(node,inode,icell)
		res .= zeros(Tv)
		src .= zeros(Tv)
		@views UK .= U[:,node.index]
		reactionwrap(res,UK)
		sourcewrap(src)

		# Loop through species
		for ispec in 1: nspecies

			# If species is enabled in the current node, add
			# integral term
			if system.node_dof[ispec,node.index]==ispec
				nf_tf = system.cellnodefactors[inode,icell]*
					tf[node.index]

				# Only add reaction term if impact ionization is
				# disabled
				if !device.model.impact
					integral[ispec] += nf_tf*(res[ispec]-
						src[ispec])
				else
					integral[ispec] += -nf_tf*src[ispec]
				end

			# if system.node_dof[ispec,node.index]==ispec
			end

		# for ispec in 1: nspecies
		end

	# for inode in 1: num_nodes(geom)
	end

	# Add contribution of flux term if field-dependent mobility is not
	# enabled
	if !device.model.fldmob
	for iedge in 1: num_edges(geom)
		VoronoiFVM._fill!(edge,iedge,icell)
		@views UKL[1: nspecies] .= U[:,edge.node[1]]
		@views UKL[nspecies+1: 2*nspecies] .= U[:,edge.node[2]]
		res .= zero(Tv)
		fluxwrap(res,UKL)
		for ispec=1: nspecies
			if system.node_dof[ispec,edge.node[1]]==ispec && 
				system.node_dof[ispec,edge.node[2]]==ispec
				integral[ispec] += system.celledgefactors[iedge,
					icell]*res[ispec]*(tf[edge.node[1]]-
					tf[edge.node[2]])
			end
		end
	end

	# if !device.model.fldmob
	end

	# for icell in 1: num_cells(grid)
	end

	return (integral)

end

# Create non-equilibrium system
# Keyword arguments:
#	unknown_storage	Sparse or dense storage for solution array
#	qfp		Use quasi-Fermi potentials for discretization
#			Current integration will not work if tf_conc is true
#	auto_sparsity	Auto-detect sparsity pattern of generic operator. This
#			is slow and should only be used for debugging
function non_equilib_sys(device,n_contacts;unknown_storage=:sparse,qfp=false,
	auto_sparsity=false)

	# Do not allow concentration test function if using QFPs
	@assert !(qfp&&tf_conc)

	# Extract model and grid from Device object
	m = device.model
	grid = device.grid
	n_boundary = length(device.b_types)

	# Define physics functions for either QFP or Boltzmann discretization
	if qfp
		flux!,source!,reaction!,bcond! = physics_noneq_qfp(device)
	else
		flux!,source!,reaction!,bcond! = physics_noneq_boltz(device)
	end

	# Pass bias and other parameters to physics functions here
	generic_noneq!(f,u,sys) = _generic_noneq!(f,u,sys,device)
	function generic_noneq_sparsity(sys)
		if !auto_sparsity
			return (_generic_noneq_sparsity(sys,device))
		else
			return (_auto_sparsity(sys,generic_noneq!))
		end
	end

	# Only use generic operator if necessary since it slows down convergence
	if m.impact && m.fldmob
		sys = VoronoiFVM.System(
			grid,
			source=source!,
			generic=generic_noneq!,
			generic_sparsity=generic_noneq_sparsity,
			breaction=bcond!,
			species=[1],
			nparams=n_contacts,
			unknown_storage=unknown_storage
		)
	elseif m.impact && !m.fldmob
		sys = VoronoiFVM.System(
			grid,
			flux=flux!,
			source=source!,
			generic=generic_noneq!,
			generic_sparsity=generic_noneq_sparsity,
			breaction=bcond!,
			species=[1],
			nparams=n_contacts,
			unknown_storage=unknown_storage
		)
	elseif !m.impact && m.fldmob
		sys = VoronoiFVM.System(
			grid,
			source=source!,
			reaction=reaction!,
			generic=generic_noneq!,
			generic_sparsity=generic_noneq_sparsity,
			breaction=bcond!,
			species=[1],
			nparams=n_contacts,
			unknown_storage=unknown_storage
		)
	else
		sys = VoronoiFVM.System(
			grid,
			flux=flux!,
			source=source!,
			reaction=reaction!,
			breaction=bcond!,
			species=[1],
			nparams=n_contacts,
			unknown_storage=unknown_storage
		)
	end

	# Carrier concentrations or QFPs only solved in semiconductor regions
	semi_regions = sort([k for (k,v) in device.r_types if 
		v=="semiconductor"])
	enable_species!(sys,2,semi_regions)
	enable_species!(sys,3,semi_regions)

	return (sys)

# function non_equilib_sys
end

# Use VFVM integration method if generic operator is not needed
function do_integrate(sys,tf,sol,device,idx)
	m = device.model
	if m.fldmob || m.impact
		return (integrate_generic(sys,tf,sol,device,idx))
	else
		return (VoronoiFVM.integrate_stdy(sys,tf,sol))
	end
end

# Solve and integrate non-equilibrium system
# Keyword arguments:
#	Plotter		specifies plotter to use for scalar plot of potential
#	verbose		verbose output for VoronoiFVM.solve!() (currently does
#	damp_initial	Initial Newton damping for solve
#	int_contacts	List of contacts to integrate
#	tf_conc		Use concentration-weighted test function for current
#			integration
#	max_round	Maximum number of iterations within roundoff error
#	max_iterations	Maximum Newton iterations
#	catch_conv	Catch ConvergenceError and break the bias loop
function non_equilib(device,bias_list,ic;Plotter=nothing,verbose=false,
	damp_initial=1.0,int_contacts=nothing,tf_conc=false,max_round=1000,
	max_iterations=100,catch_conv=false,tol_absolute=1e-10,
	tol_relative=1e-10,return_tfs=false,return_matrices=false)

	# Extract model and grid from Device object
	m = device.model
	grid = device.grid

	# 'bias_list' has as many columns as there are contacts
	n_contacts = size(bias_list,2)
	bias = zeros(n_contacts)

	# Create system
	sys = non_equilib_sys(device,n_contacts)

	# Create a solution array
	inival = unknowns(sys)
	inival .= ic
	solution = unknowns(sys)
	solution .= 0.0

	# Setup for generic operator
	idx = unknown_indices(solution)

	# Setup for flux integration
	if int_contacts!=nothing
		contact_regions = sort([k for (k,v) in device.b_types if 
			v=="contact"])
		n_contacts = length(contact_regions)
		n_integrate = length(int_contacts)
		j_list = zeros(size(bias_list,1),length(int_contacts))
	end

	# If using regular test function, define outside of bias loop to save
	# time
	if !tf_conc && (int_contacts!=nothing)
		tfs = Vector{Vector{Float64}}(undef,n_integrate)
		factory = VoronoiFVM.TestFunctionFactory(sys)
		for i in 1: n_integrate
			bc0 = [ireg for ireg in contact_regions if 
				ireg!=int_contacts[i]]
			bc1 = [int_contacts[i]]
			tfs[i] = testfunction(factory,bc0,bc1)
		end
	end

	# Define these vectors here so we can return them
	if tf_conc && (int_contacts!=nothing)
		tfs_n = Vector{Vector{Float64}}(undef,n_integrate)
		tfs_p = Vector{Vector{Float64}}(undef,n_integrate)
	end

	# Matrix accumulator
	matrix_list = []

	# Loop over bias voltages here to avoid repeating discretization
	for i in 1: size(bias_list,1)
		bias = bias_list[i,: ]

		# Using max_round here can help convergence errors when the IC
		# is very close to the solution, but should not be used in
		# general
		print(bias," ")
		control = VoronoiFVM.SolverControl()
		control.verbose = verbose
		control.damp_initial = damp_initial
		control.max_round = max_round
		control.max_iterations = max_iterations
		control.tol_absolute=tol_absolute
		control.tol_relative=tol_relative
		if !catch_conv
			solve!(solution,inival,sys,control=control,params=bias)
		else
			try
				solve!(solution,inival,sys,control=control,
					params=bias)
			catch e
				if isa(e,VoronoiFVM.ConvergenceError)
					println("Caught ConvergenceError")
					break
				else
					rethrow(e)
				end
			end
		end
		inival .= solution

		# Compute condition number if requested and add it to
		# accumulator
		if return_matrices

			# Remove penalty factors from the Jacobian
			m = Matrix(sys.matrix)
			ident = one(zeros(3,3))
			m[1: 3,: ] .= 0.0
			m[end-2: end,: ] .= 0.0
			m[1: 3,1: 3] = ident
			m[end-2: end,end-2: end] = ident

			# Push to accumulator
			push!(matrix_list,m)

		end

		# Use concentration-weighted test function
		# This is slow since the test function must be recalculated for
		# each bias point
		if tf_conc && (int_contacts!=nothing)

			# Create vectors to hold test functions (2 for each
			# contact, for electron and hole currents) and their
			# initial conditions
			ic_tf_n = zeros(Float64,1,size(solution,2))
			ic_tf_p = zeros(Float64,1,size(solution,2))

			# Loop through contacts and create test functions
			for i in 1: n_integrate
				bc0 = [ireg for ireg in contact_regions if 
					ireg!=int_contacts[i]]
				bc1 = [int_contacts[i]]
				factory_n = tf_sys(sys,solution,device,bc0,bc1,
					jp=false)
				factory_p = tf_sys(sys,solution,device,bc0,bc1,
					jp=true)
				tfs_n[i] = tf_solve(factory_n,ic_tf_n)
				tfs_p[i] = tf_solve(factory_p,ic_tf_p)
				ic_tf_n .= tfs_n[i]'
				ic_tf_p .= tfs_p[i]'
			end

		end

		# Plot solution if Plotter is specified in args
		scalarplot(grid,solution[1,: ],title="Non-equilibrium \
			potential",Plotter=Plotter,colormap=:coolwarm,
			resolution=(1000,500))

		# Get total current at each contact
		if int_contacts!=nothing for j in 1: n_integrate
			print(j," ")
			if tf_conc
				bflux_n = do_integrate(
					sys,tfs_n[j],solution,device,idx)
				bflux_p = do_integrate(
					sys,tfs_p[j],solution,device,idx)
				j_list[i,j] = m.q*(bflux_n[3]-bflux_p[2])
			else
				bflux = do_integrate(
					sys,tfs[j],solution,device,idx)
				j_list[i,j] = m.q*(bflux[3]-bflux[2])
			end
		end end

	end
	println()

	# Return value logic
	r = [nothing,sys,solution]
	if int_contacts!=nothing
		r[1] = j_list
	end
	if return_tfs
		if tf_conc
			push!(r,tfs_n)
			push!(r,tfs_p)
		else
			push!(r,tfs)
		end
	end
	if return_matrices
		push!(r,matrix_list)
	end
	return (tuple(r...))

# function non_equilib
end

