# I-V curve generation via arc length continuation
# Sam Chinnery
# 2022-02-05

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

# Contains Newton state and solver parameters
mutable struct IVState{Tz,Tv,Ti}
	z::AbstractVector{Tz}		# New solution (updated inplace)
	v::Tv				# New voltage (updated inplace)
	i::Tv				# New current (updated inplace)
	zpred::AbstractVector{Tz}	# Predicted solution
	vpred::Tv			# Predicted voltage
	ipred::Tv			# Predicted current
	z0::AbstractVector{Tz}		# Old solution
	v0::Tv				# Old voltage
	i0::Tv				# Old current
	N::Tv				# Arc length residual
	bias::AbstractVector{Tv}	# Bias vector
	c_active::Ti			# Active contact
	dfdv::AbstractVector{Tv}	# Vector derivative dF/dV
	didz::AbstractVector{Tv}	# Gradient vector dI/dz
	dids::Tv			# Scalar derivative dI/ds
	dvds::Tv			# Scalar derivative dV/ds
	dzds::AbstractVector{Tv}	# Vector derivative dz/ds
	dndz::AbstractVector{Tv}	# Gradient vector dN/dz
	dndv::Tv			# Scalar derivative dN/dV
	theta::Tv			# Weight parameter
	ds::Tv				# Arc length
	use_n1::Bool			# Residual function choice
	log_j::Bool			# Use log10 of current
	scale_j::Tv			# Scale factor for current
	tol_abs::Tv			# Absolute Newton tolerance
	tol_rel::Tv			# Relative Newton tolerance
	tol_mono::Tv			# Monotonicity Newton tolerance
	iteration::Ti			# Current iteration number
	damp_initial::Tv		# Initial damping factor
	damp::Tv			# Newton damping factor
	n_damp::Ti			# Number of damped iterations
	max_iters::Ti			# Maximum Newton iterations
	du::Tv				# Current update norm
	du_prev::Tv			# Previous update norm
	du_first::Tv			# First update norm
	converged::Ti			# Newton convergence status
	verbose::Ti			# Print debug information
end

# Contains IV curve and continuation parameters
mutable struct IVCurve{Tv,Ti}
	v::AbstractVector{Tv}		# Holds voltage
	i::AbstractVector{Tv}		# Holds current
	vmin::Tv			# Minimum bias voltage
	vmax::Tv			# Maximum bias voltage
	tol_trunc::Tv			# Truncation error tolerance
	max_steps::Ti			# Maximum continuation steps
	step::Ti			# Current step number
	ds_decrease_factor::Tv		# Factor to decrease ds by
	ds_max_increase::Tv		# Maximum increase of ds
	trunc_error_norm::Tv		# Norm of truncation error
	dsmin::Tv			# Minimum ds
	n_solve::Ti			# Number of Newton solves
	n_solve_max::Ti			# Maximum Newton solves
	skip_restart::Bool		# Whether to skip restart
end

# Default constructor
function IVState(
	ic::AbstractVector{Tz},
	bias::AbstractVector{Tv},
	c_active::Ti,
	ds::Tv
) where {Tz,Tv,Ti}

	dval = Tv(0.0)
	dint = Ti(0)

	v = bias[c_active]
	s = IVState{Tz,Tv,Ti}(
		copy(ic),copy(v),copy(dval),
		copy(ic),copy(v),copy(dval),
		copy(ic),copy(v),copy(dval),
		dval,copy(bias),c_active,
		zeros(Tv,0),zeros(Tv,0),dval,dval,zeros(Tv,0),zeros(Tv,0),dval,
		Tv(1.0),ds,true,false,Tv(1.0),Tv(1e-10),Tv(1e-10),Tv(1e-3),
		Ti(1),Tv(1.0),Tv(1.0),Ti(0),Ti(10),dval,dval,dval,
		Ti(status_not_converged),Ti(verbose_noverbose)
	)
	return (s)

# function IVState
end

# Default constructor
function IVCurve(
	state::IVState{Tz,Tv,Ti},
	vmin::Tv,
	vmax::Tv,
	dsmin::Tv
) where {Tz,Tv,Ti}

	s = IVCurve{Tv,Ti}(
		zeros(Tv,0),zeros(Tv,0),vmin,vmax,
		Tv(0.1),Ti(10),Ti(1),Tv(0.5),Tv(2.0),Tv(0.0),
		dsmin,Tv(0.0),Ti(20),false
	)
	return (s)

# function IVCurve
end

# Initialize Newton state
function initialize!(
	state::IVState{Tz,Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	state.converged = status_not_converged
	state.iteration = 1
	state.damp = state.damp_initial

	state.du = 0
	state.du_prev = 0
	state.du_first = 0

	state.z0 = state.z
	state.v0 = state.v

	current!(state,d,sys,tf)
	state.i0 = state.i
	update_res_jac!(state,sys)
	assemble_derivs!(state,d,sys,tf,tangent=true)

	return (nothing)

# function initialize!
end

# Compute derivative of residual with respect to active contact voltage
function assemble_dfdv(
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv,Ti2},
	c_active::Ti
) where {Tv,Ti,Ti2}

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
			if reg==c_active
				bnode_idx = idx[1,bnode.index]
				dfdv[bnode_idx] = -VoronoiFVM.Dirichlet
			end
		end

	# for ibnode in 1: bnodes_per_bface
	end

	return (dfdv)

# function assemble_dfdv
end

# Inplace dF/dV
function assemble_dfdv!(
	state::IVState{Tz,Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv,Ti2}
) where {Tz,Tv,Ti,Ti2}

	dfdv = assemble_dfdv(d,sys,state.c_active)
	state.dfdv = dfdv
	return (nothing)

# function assemble_dfdv!
end

# Gradient of current with respect to solution
function assemble_didz!(
	state::IVState{Tz,Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	i_z(z0) = current(z0,state,d,sys,tf)
	state.didz = ReverseDiff.gradient(i_z,state.z)
	return (nothing)

# function assemble_didz!
end

# Gradient of N with respect to solution, and derivative of N with respect to V
function assemble_dndz_dndv!(
	state::IVState{Tz,Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	N_z(z) = arc_constraint(z,state.v,state,d,sys,tf)
	N_v(v) = arc_constraint(state.z,v,state,d,sys,tf)
	state.dndz = ReverseDiff.gradient(N_z,state.z)
	state.dndv = ForwardDiff.derivative(N_v,state.v)
	return (nothing)

# function assemble_dndz_dndv!
end

# Assemble all derivatives into state object
function assemble_derivs!(
	state::IVState{Tz,Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv};
	tangent=false
) where {Tz,Tv,Ti}

	# dF/dV and dI/dz do not depend on any other derivatives
	assemble_dfdv!(state,d,sys)
	assemble_didz!(state,d,sys,tf)

	# dV/ds, dz/ds and dI/ds depend on dF/dv and dI/dz
	if tangent
		get_tangent!(state,d,sys,tf)
	end

	# N, dN/dz and dN/dV depend on dV/ds and dI/ds
	arc_constraint!(state,d,sys,tf)
	assemble_dndz_dndv!(state,d,sys,tf)

	# Print debugging information if requested
	print_derivs(state)

	return (nothing)

# function assemble_derivs!
end

# Integrate using test function to find current
# This is used by other functions that do not pass state::IVState
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

# Integrate using test function to find current
# This does not update the solution in the state object
function current(
	z,
	state::IVState{Tz,Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	j = state.scale_j*current(z,d,sys,tf)
	if state.log_j
		return (log10(abs(j)))
	else
		return (j)
	end

# function current
end

# Integrate current and update the state object
function current!(
	state::IVState{Tz,Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	j = current(state.z,state,d,sys,tf)
	state.i = j
	return (nothing)

# function current!
end

# Computes the tangent vector without updating state, for AD purposes
function get_tangent!(
	state::IVState{Tz,Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	jinv_dfdv = lu(-sparse(sys.matrix))\state.dfdv

	dvds = (1+dot(state.didz,jinv_dfdv)^2)^-0.5
	dzds = jinv_dfdv*dvds
	dids = dot(state.didz,dzds)

	if state.iteration==1
		sign_dvds = 1.0
	else
		sign_dvds = sign(dot([dids,dvds],[state.dids,state.dvds]))
	end
	state.dvds = sign_dvds*dvds
	state.dzds = sign_dvds*dzds
	state.dids = sign_dvds*dids

	return (nothing)

# function get_tangent
end

# Compute arc length constraint N(I(s),V(s))
# This does not update N in the state object
function arc_constraint(
	z,
	v,
	state::IVState{Tz,Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	# Using non-inplace calls to maintain differentiability
	j = current(z,state,d,sys,tf)

	# Compute either N1 or N2 from Coughran et al. 1989
	delta_i = j-state.i0
	delta_v = v-state.v0
	if state.use_n1
		N = state.theta*state.dids*delta_i 
		N += (2-state.theta)*state.dvds*delta_v - state.ds
	else
		N = delta_i^2 + delta_v^2 - state.ds^2
	end
	return (N)

# function arc_constraint
end

# Compute arc length constraint and update state object
function arc_constraint!(
	state::IVState{Tz,Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	N = arc_constraint(state.z,state.v,state,d,sys,tf)
	state.N = N
	return (nothing)

# function arc_constraint!
end

# Update residual and Jacobian inplace
function update_res_jac!(
	state::IVState{Tz,Tv,Ti},
	sys::VoronoiFVM.AbstractSystem{Tv}
) where {Tz,Tv,Ti}

	uhash = hash((state.z,state.bias))
	if uhash!=sys.uhash
		state.bias[state.c_active] = state.v
		ur = reshape(state.z,sys)
		VoronoiFVM.eval_and_assemble(
			sys,ur,ur,sys.residual,Inf,Inf,0.0,state.bias)
		sys.uhash = uhash
		sys.history.nd += 1
	end
	return (nothing)

# function update_res_jac!
end

# Print derivative information for debugging
function print_derivs(
	state::IVState
)

	if state.verbose < verbose_print_derivs
		return (nothing)
	end

	println("Derivative assembly")

	nnz = count(v->v!=0.0,state.dfdv)
	dnorm = norm(state.dfdv)
	println("\tdF/dV: ",dnorm," norm, ",nnz," nonzero")

	nnz = count(v->v!=0.0,state.didz)
	dnorm = norm(state.didz)
	println("\tdI/dz: ",dnorm," norm, ",nnz," nonzero")

	nnz = count(v->v!=0.0,state.dzds)
	dnorm = norm(state.dzds)
	println("\tdV/ds = ",state.dvds)
	println("\tdz/ds: ",dnorm," norm, ",nnz," nonzero")
	println("\tdI/ds = ",state.dids)

	nnz = count(v->v!=0.0,state.dndz)
	dnorm = norm(state.dndz)
	println("\tdN/dz: ",dnorm," norm, ",nnz," nonzero")
	println("\tdN/dV = ",state.dndv)

	return (nothing)

# function print_derivs
end

# Print Newton information for debugging
function print_newton(
	state::IVState,
	sys::VoronoiFVM.AbstractSystem{Tv},
	update::AbstractVector{Tv}
) where {Tv}

	if state.verbose < verbose_print_iters
		return (nothing)
	end

	nnz_F = count(v->v!=0.0,sys.residual)
	norm_F = norm(sys.residual)
	println("Iteration ",state.iteration)
	println("\tF: ",norm_F," norm, ",nnz_F," nonzero")
	println("\tN = ",state.N)

	nnz_update = count(v->v!=0.0,update)
	norm_update = norm(update)
	max_update = norm(update,Inf)
	print("\tupdate: ",norm_update," norm, ",nnz_update," nonzero, ")
	println(max_update," max")

	println("\tI = ",state.i)
	println("\tV = ",state.v)

# function print_newton
end

# Predict next solution from tangent vector
function predict!(
	state::IVState
)

	state.z += state.ds*state.dzds
	#state.v += state.ds*state.dvds
	state.i += state.ds*state.dids
	state.zpred = state.z
	state.vpred = state.v
	state.ipred = state.i
	return (nothing)

# function predict!
end

# Check if Newton solver has converged
function check_convergence!(
	state::IVState{Tz,Tv,Ti}
) where {Tz,Tv,Ti}

	if state.iteration==1
		state.du_first = state.du
		state.du_prev = state.du
	end

	# Absolute or relative tolerance met
	if state.du < state.tol_abs
		state.converged = status_converged
		return (nothing)
	end
	if state.du/state.du_first < state.tol_rel
		state.converged = status_converged
		return (nothing)
	end

	# Monotonicity error
	if state.du/state.du_prev > 1/state.tol_mono
		state.converged = err_monotonicity
		return (nothing)
	end

	# Max iterations reached
	if state.iteration >= state.max_iters
		state.converged = err_max_iters
		return (nothing)
	else
		state.iteration += 1
	end

	state.du_prev = state.du
	return (nothing)

# function check_convergence!
end

# Compute new damping value
function update_damp!(
	state::IVState{Tz,Tv,Ti}
) where {Tz,Tv,Ti}

	if state.iteration > state.n_damp
		state.damp = 1.0
	end
	return (nothing)

# function update_damp!
end

# Single damped Newton step
function newton_step!(
	state::IVState{Tz,Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	J = sparse(sys.matrix)
	F = VoronoiFVM.values(sys.residual)
	A = vcat(hcat(J,state.dfdv),hcat(state.dndz',state.dndv))
	b = vcat(F,state.N)
	update = lu(-A)\b
	update_z = update[1: end-1]
	update_v = update[end]

	update_damp!(state)
	state.z += state.damp*update_z
	state.v += state.damp*update_v
	current!(state,d,sys,tf)
	state.du = norm(update,Inf)

	print_newton(state,sys,update)
	return (nothing)

# function newton_step!
end

# Newton process
function newton!(
	state::IVState{Tz,Tv,Ti},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	while state.converged==status_not_converged
		update_res_jac!(state,sys)
		assemble_derivs!(state,d,sys,tf)
		newton_step!(state,d,sys,tf)
		check_convergence!(state)
	end
	return (nothing)

# function newton!
end

# Print debug information for continuation step
function print_cont_1(
	state::IVState{Tz,Tv,Ti},
	curve::IVCurve{Tv}
) where {Tz,Tv,Ti}

	if state.verbose < verbose_print_cont
		return (nothing)
	end

	printstyled("Continuation step ",curve.step,color=:light_yellow)
	println(", ds = ",state.ds)
	println("\tVinitial = ",state.v0,", Iinitial = ",state.i0)
	println("\tVpred = ",state.vpred,", Ipred = ",state.ipred)
	println("\tVfinal = ",state.v,", Ifinal = ",state.i)
	println("\tNewton iterations: ",state.iteration)

	if state.converged==status_converged
		convergence_color = :green
	else
		convergence_color = :light_red
	end
	printstyled("\tConvergence status: ",color=convergence_color)
	printstyled(state.converged,"\n",color=convergence_color)
	return (nothing)

# function print_cont_1
end

# Print debug information after continuation update
function print_cont_2(
	state::IVState{Tz,Tv,Ti},
	curve::IVCurve{Tv}
) where {Tz,Tv,Ti}

	if state.verbose < verbose_print_cont
		return (nothing)
	end

	if curve.trunc_error_norm>curve.tol_trunc
		trunc_color = :light_red
	else
		trunc_color = :green
	end
	printstyled("\tTruncation error: ",color=trunc_color)
	printstyled(curve.trunc_error_norm,"\n\n",color=trunc_color)
	return (nothing)

# function print_cont_2
end

# Reset state to the previous continuation step
function reset_state!(
	state::IVState{Tz,Tv,Ti}
) where {Tz,Tv,Ti}

	state.z = state.z0
	state.v = state.v0
	state.i = state.i0
	return (nothing)

# function reset_state!
end

# Update step size
function ds_update!(
	state::IVState{Tz,Tv,Ti},
	curve::IVCurve{Tv},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	trunc_error = [state.i-state.ipred,state.v-state.vpred]
	curve.trunc_error_norm = norm(trunc_error)

	if state.converged==status_converged
		ds_factor = sqrt(curve.tol_trunc/curve.trunc_error_norm)
		ds_factor_clamp = min(ds_factor,curve.ds_max_increase)
		if curve.trunc_error_norm>curve.tol_trunc
			if curve.skip_restart
				curve.skip_restart = false
				curve.step += 1
			else
				state.ds *= ds_factor_clamp
				reset_state!(state)
				curve.skip_restart = true
			end
		else
			state.ds *= ds_factor_clamp
			curve.step += 1
		end
	else
		reset_state!(state)
		state.ds *= curve.ds_decrease_factor
	end
	return (nothing)

# function ds_update!
end

# Check if continuation should be terminated
function check_cont_terminate(
	state::IVState{Tz,Tv,Ti},
	curve::IVCurve{Tv}
) where {Tz,Tv,Ti}

	if curve.n_solve >= curve.n_solve_max
		error("max solves reached")
	end
	if abs(state.ds) < curve.dsmin
		error("dsmin reached")
	end

	curve.n_solve += 1
	return (nothing)

# function check_cont_terminate
end

# Continuation
function continuation!(
	state::IVState{Tz,Tv,Ti},
	curve::IVCurve{Tv},
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem{Tv},
	tf::AbstractArray{Tv}
) where {Tz,Tv,Ti}

	while curve.step<=curve.max_steps
		check_cont_terminate(state,curve)
		initialize!(state,d,sys,tf)
		predict!(state)
		newton!(state,d,sys,tf)
		print_cont_1(state,curve)
		ds_update!(state,curve,d,sys,tf)
		print_cont_2(state,curve)
		push!(curve.v,state.v)
		push!(curve.i,state.i)
	end
	return (nothing)

# function continuation!
end

