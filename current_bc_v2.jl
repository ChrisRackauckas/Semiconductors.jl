# Current boundary conditions
# Sam Chinnery
# 2022-02-16

# Verbosity values
const verbose_noverbose = 0
const verbose_print_cont = 1
const verbose_print_iters = 2
const verbose_print_derivs = 3
const verbose_print_all = typemax(Int)

# Convergence errors
const status_converged = 1
const status_not_converged = 0
const err_max_iters = -1
const err_monotonicity = -2

# Contains Newton state
mutable struct NewtonState{Tz,Tv,Ti}

	z::AbstractVector{Tz}		# New solution (updated inplace)
	bias::AbstractVector{Tv}	# Bias vector (updated inplace)
	i::Tv				# New current (updated inplace)
	dfdv::AbstractVector{Tv}	# Vector derivative dF/dV
	didz::AbstractVector{Tv}	# Gradient vector dI/dz
	didv::Tv			# Scalar derivative dI/dV (usually 0.0)
	iteration::Ti			# Current iteration number
	damp::Tv			# Newton damping factor
	update_z::AbstractVector{Tv}	# Update to solution
	update_v::Tv			# Update to bias voltage
	du::Tv				# Current update norm
	du_prev::Tv			# Previous update norm
	du_first::Tv			# First update norm
	converged::Ti			# Newton convergence status

# mutable struct NewtonState
end

# Contains Newton parameters
mutable struct NewtonParams{Tv,Ti}

	set_current::Tv			# Current BC value
	c_active::Ti			# Active contact
	tol_abs::Tv			# Absolute Newton tolerance
	tol_rel::Tv			# Relative Newton tolerance
	tol_mono::Tv			# Monotonicity Newton tolerance
	max_iters::Ti			# Maximum Newton iterations
	do_damp_search::Bool		# Search for initial damping value
	damp_initial::Tv		# Default initial damping
	damp_growth::Tv			# Damping update parameter
	damp_search_iters::Ti		# Maximum iterations for damping search
	damp_search_decrease::Tv	# Damping decrease for search
	dirichlet_scale::Tv		# Scale factor for BC assembly
	verbose::Ti			# Print debug information

# mutable struct NewtonParams
end

# Default constructor
function NewtonState(
	ic::AbstractVector{Tz},
	bias::AbstractVector{Tv}
) where {Tz,Tv}

	s = NewtonState(
		ic,bias,Tv(0.0),zeros(Tv,0),zeros(Tv,0),Tv(0.0),1,
		Tv(1.0),zeros(Tv,0),Tv(0.0),Tv(0.0),Tv(0.0),Tv(0.0),
		status_not_converged
	)
	return (s)

# function NewtonState
end

# Default constructor
function NewtonParams(
	set_current::Tv,
	c_active::Ti;
	tol_abs=1e-10,
	tol_rel=1e-10,
	tol_mono=1e-3,
	max_iters=25,
	do_damp_search=false,
	damp_initial=0.1,
	damp_growth=1.5,
	damp_search_iters=10,
	damp_search_decrease=0.5,
	dirichlet_scale=1.0,
	verbose=verbose_noverbose
) where {Tv,Ti}

	p = NewtonParams{Tv,Ti}(
		set_current,c_active,tol_abs,tol_rel,tol_mono,max_iters,
		do_damp_search,damp_initial,damp_growth,damp_search_iters,
		damp_search_decrease,dirichlet_scale,verbose
	)
	return (p)

# function NewtonParams
end

# Initialize Newton state
function initialize!(
	state::NewtonState{Tz,Tv,Ti},
	params::NewtonParams{Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv}
) where {Tz,Tv,Ti}

	state.converged = status_not_converged
	state.iteration = 1
	state.damp = params.damp_initial

	state.du = 0
	state.du_prev = 0
	state.du_first = 0

# function initialize!
end

# Integrate using test function to find current
# This is used by other functions that do not pass state::NewtonState
function current(
	z,
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tv}

	# Perform integration
	idx = unknown_indices(unknowns(sys))
	zr = reshape(z,sys)
	bflux = do_integrate(sys,tf,zr,d,idx)
	j = d.model.q*(bflux[3]-bflux[2])
	return (j)

# function current
end

# Integrate using test function and update the NewtonState
function current!(
	state::NewtonState{Tz,Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	state.i = current(state.z,d,sys,tf)
	return (nothing)

# function current!
end

# Update residual and Jacobian
function update_res_jac!(
	state::NewtonState{Tz,Tv,Ti},
	sys::VoronoiFVM.AbstractSystem{Tv}
) where {Tz,Tv,Ti}

	uhash = hash((state.z,state.bias))
	if uhash!=sys.uhash
		zr = reshape(state.z,sys)
		VoronoiFVM.eval_and_assemble(
			sys,zr,zr,sys.residual,Inf,Inf,0.0,state.bias)
		sys.uhash = uhash
		sys.history.nd += 1
	end
	return (nothing)

# function update_res_jac!
end

# Compute derivative of residual with respect to active contact voltage
function assemble_dfdv!(
	state::NewtonState{Tz,Tv,Ti},
	params::NewtonParams{Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv,Ti2}
) where {Tz,Tv,Ti,Ti2}

	# Get geometry information from grid
	n_bfaces = num_bfaces(d.grid)
	bgeom = d.grid[BFaceGeometries][1]
	bnodes_per_bface = num_nodes(bgeom)
	bnode = VoronoiFVM.BNode{Tv,Ti2}(sys)

	# Create vector to hold dF/dv and get unknown indices within solution
	# vector
	dfdv = zeros(Tv,num_dof(sys))
	idx = unknown_indices(unknowns(sys))

	# Loop through all boundary faces
	for ibface in 1: n_bfaces

		# Get boundary region of current boundary face
		reg = d.grid[BFaceRegions][ibface]

		# Loop through all boundary nodes in each boundary face and add
		# Dirichlet penalty factor if bnode is in the active contact
		for ibnode in 1: bnodes_per_bface
			VoronoiFVM._fill!(bnode,ibnode,ibface)
			if reg==params.c_active
				bnode_idx = idx[1,bnode.index]
				dfdv[bnode_idx] = -VoronoiFVM.Dirichlet
			end
		end

	# for ibnode in 1: bnodes_per_bface
	end

	state.dfdv = dfdv
	return (nothing)

# function assemble_dfdv!
end

# Gradient of current with respect to solution
function assemble_didz!(
	state::NewtonState{Tz,Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	i_z(z0) = current(z0,d,sys,tf)
	state.didz = ReverseDiff.gradient(i_z,state.z)
	return (nothing)

# function assemble_didz!
end

# Check if Newton is converged
function check_convergence!(
	state::NewtonState{Tz,Tv,Ti},
	params::NewtonParams{Tv,Ti}
) where {Tz,Tv,Ti}

	if state.iteration==1
		state.du_first = state.du
		state.du_prev = state.du
	end

	# Absolute or relative tolerance met
	if state.du < params.tol_abs
		state.converged = status_converged
		return (nothing)
	end
	if state.du/state.du_first < params.tol_rel
		state.converged = status_converged
		return (nothing)
	end

	# Monotonicity error
	if state.du/state.du_prev > 1/params.tol_mono
		state.converged = err_monotonicity
		error("Monotonicity error")
	end

	# Max iterations reached
	if state.iteration >= params.max_iters
		state.converged = err_max_iters
		error("Max iterations reached")
	else
		state.iteration += 1
	end

	state.du_prev = state.du
	state.converged = status_not_converged
	return (nothing)

# function check_convergence
end

# Assemble derivatives and get current
function assemble_derivs!(
	state::NewtonState{Tz,Tv,Ti},
	params::NewtonParams{Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	update_res_jac!(state,sys)
	assemble_dfdv!(state,params,d,sys)
	assemble_didz!(state,d,sys,tf)
	current!(state,d,sys,tf)

# function assemble_derivs!
end

# Print Newton information for debugging
function print_newton_current(
	state::NewtonState{Tz,Tv,Ti},
	params::NewtonParams{Tv,Ti},
	sys::VoronoiFVM.AbstractSystem{Tv},
	update::AbstractVector{Tv}
) where {Tz,Tv,Ti}

	if params.verbose < verbose_print_iters
		return (nothing)
	end

	nnz_F = count(v->v!=0.0,sys.residual)
	norm_F = norm(sys.residual)
	res_current = state.i/params.set_current-1.0
	println("Iteration ",state.iteration)
	println("\tF: ",norm_F," norm, ",nnz_F," nonzero")
	println("\tres_current = ",res_current)

	nnz_update = count(v->v!=0.0,update)
	norm_update = norm(update)
	max_update = norm(update,Inf)
	print("\tupdate: ",norm_update," norm, ",nnz_update," nonzero, ")
	println(max_update," max")
	println("\tdamp = ",state.damp)

	println("\tI = ",state.i)
	println("\tbias: ",state.bias)

	return (nothing)

# function print_newton_current
end

# Search for initial damping if default value is too high
function damp_search!(
	state::NewtonState{Tz,Tv,Ti},
	params::NewtonParams{Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	# Compute initial residual
	F = VoronoiFVM.values(sys.residual)
	res_current = state.i/params.set_current-1.0
	b = vcat(F,res_current)
	res_norm_old = norm(b,Inf)
	z_old = copy(state.z)
	bias_old = copy(state.bias)
	damp_old = state.damp

	# Search for a fixed number of iterations to avoid excessively small
	# damping
	for j in 1: params.damp_search_iters

		# Evaluate the residual with the candidate damping value
		state.z -= state.damp*state.update_z
		state.bias[params.c_active] -= state.damp*state.update_v
		update_res_jac!(state,sys)
		current!(state,d,sys,tf)
		F = VoronoiFVM.values(sys.residual)
		res_current = state.i/params.set_current-1.0
		b = vcat(F,res_current)
		res_norm_new = norm(b,Inf)

		# If the candidate damping does not sufficiently decrease the
		# residual, decrease it and try again
		if res_norm_new > res_norm_old
			state.damp *= params.damp_search_decrease
		else
			break
		end
		state.z .= z_old
		state.bias .= bias_old

	# for j in 1: params.damp_search_iters
	end

	# Reset state variables, since newton_step!() will update them anyway
	state.z .= z_old
	state.bias .= bias_old
	return (nothing)

# function damp_search!
end

# Update damping value
function update_damp!(
	state::NewtonState{Tz,Tv,Ti},
	params::NewtonParams{Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	if state.iteration==1
		if params.do_damp_search
			damp_search!(state,params,d,sys,tf)
		end
	else
		state.damp *= params.damp_growth
		state.damp = min(1.0,state.damp)
	end
	return (nothing)

# function update_damp!
end

# Single Newton step
function newton_step_current!(
	state::NewtonState{Tz,Tv,Ti},
	params::NewtonParams{Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	J = sparse(sys.matrix)
	F = VoronoiFVM.values(sys.residual)
	res_current = state.i/params.set_current-1.0
	didz_scaled = 1/params.set_current * state.didz'
	A = vcat(hcat(J,state.dfdv),hcat(didz_scaled,state.didv))
	b = vcat(F,res_current)

	update = lu(A)\b
	state.update_z = update[1: end-1]
	state.update_v = update[end]
	state.du = norm(update,Inf)

	update_damp!(state,params,d,sys,tf)
	state.z -= state.damp*state.update_z
	state.bias[params.c_active] -= state.damp*state.update_v

	print_newton_current(state,params,sys,update)
	return (nothing)

# function newton_step!
end

# Newton solver with current BC
function newton_current!(
	state::NewtonState{Tz,Tv,Ti},
	params::NewtonParams{Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	initialize!(state,params,d,sys)
	while state.converged==status_not_converged
		assemble_derivs!(state,params,d,sys,tf)
		newton_step_current!(state,params,d,sys,tf)
		check_convergence!(state,params)
	end

	current!(state,d,sys,tf)
	return (nothing)

# function newton_current!
end

