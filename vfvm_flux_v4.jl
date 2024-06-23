# VFVM interface to DiffEqFlux.jl
# Sam Chinnery
# 2022-02-21

# Derivative of Bernoulli function
function fbernoulli_pm_deriv(x::Real)

	dbp = ForwardDiff.derivative(fbernoulli,x)
	dbm = -1.0-dbp
	return (dbp,dbm)

# function fbernoulli_pm_deriv
end

# Replace solution with charge neutrality IC, assuming doping is defined as an
# Array instead of a function and that ni is constant
function charge_neutrality_ic!(
	state::NewtonState,
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem
)

	m = d.model
	ni = m.ni[1]
	n0p0 = m.doping_to_np.(m.doping,ni)
	n0 = [x[1] for x in n0p0]
	p0 = [x[2] for x in n0p0]
	v = m.np_to_v.(n0,p0,ni,m.vt)
	ic_mat = Matrix([v n0 p0]')
	ic = VoronoiFVM.values(ic_mat)
	state.z = ic
	return (nothing)

# function charge_neutrality_ic!(
end

# Flux part of residual and Jacobian
function assemble_flux!(
	F::Zygote.Buffer,
	J::Zygote.Buffer,
	state::NewtonState,
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem
)

	m = d.model
	idx = LinearIndices(unknowns(sys))
	dims = dim_space(d.grid)
	z = state.z

	for icell in 1: num_cells(d.grid)

		region = d.grid.components[CellRegions][icell]
		eps = m.e0*m.er[region]
		un = m.mobility_n(0,0,0.0,region,false)
		up = m.mobility_p(0,0,0.0,region,false)

		for iedge in 1: num_edges(d.grid.components[CellGeometries][1])

			# We cannot use the typical CellEdges, EdgeNodes
			# adjacencies here since they mutate the grid
			if dims==1
				cell_nodes = d.grid.components[CellNodes]
				coords = d.grid.components[Coordinates]
				inode_k = cell_nodes[1,icell]
				inode_l = cell_nodes[2,icell]
				x_k = coords[1,inode_k]
				x_l = coords[1,inode_l]
				edge_factor = 1/abs(x_l-x_k)
			elseif dims==2
				cell_faces = d.grid.components[CellFaces]
				cell_edge = cell_faces[iedge,icell]
				face_nodes = d.grid.components[FaceNodes]
				inode_k = face_nodes[1,cell_edge]
				inode_l = face_nodes[2,cell_edge]
				edge_factor = sys.celledgefactors[iedge,icell]
			end

			idx_vk = idx[1,inode_k]
			idx_vl = idx[1,inode_l]
			idx_nk = idx[2,inode_k]
			idx_nl = idx[2,inode_l]
			idx_pk = idx[3,inode_k]
			idx_pl = idx[3,inode_l]

			deltav = (z[idx_vl]-z[idx_vk])/m.vt
			bp,bm = fbernoulli_pm(deltav)
			dbp,dbm = fbernoulli_pm_deriv(deltav)

			v_flux = eps*(z[idx_vk]-z[idx_vl])
			n_flux = un*m.vt*(z[idx_nk]*bm-z[idx_nl]*bp)
			p_flux = up*m.vt*(z[idx_pk]*bp-z[idx_pl]*bm)
			v_flux *= edge_factor
			n_flux *= edge_factor
			p_flux *= edge_factor

			F[idx_vk] += v_flux
			F[idx_vl] -= v_flux
			F[idx_nk] += n_flux
			F[idx_nl] -= n_flux
			F[idx_pk] += p_flux
			F[idx_pl] -= p_flux

			j_vk_vk = edge_factor*eps

			J[idx_vk,idx_vk] += j_vk_vk
			J[idx_vk,idx_vl] -= j_vk_vk
			J[idx_vl,idx_vk] -= j_vk_vk
			J[idx_vl,idx_vl] += j_vk_vk

			j_nk_nk = edge_factor*un*m.vt*bm
			j_nl_nl = edge_factor*un*m.vt*bp

			J[idx_nk,idx_nk] += j_nk_nk
			J[idx_nk,idx_nl] -= j_nl_nl
			J[idx_nl,idx_nk] -= j_nk_nk
			J[idx_nl,idx_nl] += j_nl_nl

			j_nk_vk = un*(z[idx_nk]*dbm+z[idx_nl]*dbp)
			j_nk_vk *= edge_factor

			J[idx_nk,idx_vk] += j_nk_vk
			J[idx_nk,idx_vl] -= j_nk_vk
			J[idx_nl,idx_vk] -= j_nk_vk
			J[idx_nl,idx_vl] += j_nk_vk

			j_pk_pk = edge_factor*up*m.vt*bp
			j_pl_pl = edge_factor*up*m.vt*bm

			J[idx_pk,idx_pk] += j_pk_pk
			J[idx_pk,idx_pl] -= j_pl_pl
			J[idx_pl,idx_pk] -= j_pk_pk
			J[idx_pl,idx_pl] += j_pl_pl

			j_pk_vk = up*(z[idx_pk]*dbp+z[idx_pl]*dbm)
			j_pk_vk *= edge_factor

			J[idx_pk,idx_vk] -= j_pk_vk
			J[idx_pk,idx_vl] += j_pk_vk
			J[idx_pl,idx_vk] += j_pk_vk
			J[idx_pl,idx_vl] -= j_pk_vk

		# for iedge in 1: num_edges(d.grid[CellGeometries][1])
		end

	# for icell in 1: num_cells(d.grid)
	end

	return (nothing)

# function assemble_flux!
end

# Reaction part of residual
function assemble_reaction!(
	F::Zygote.Buffer,
	J::Zygote.Buffer,
	state::NewtonState,
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem
)

	m = d.model
	idx = LinearIndices(unknowns(sys))
	dims = dim_space(d.grid)
	z = state.z

	for icell in 1: num_cells(d.grid)

		region = d.grid.components[CellRegions][icell]

		for inode in 1: num_nodes(d.grid.components[CellGeometries][1])

			cell_nodes = d.grid.components[CellNodes]
			inode_k = cell_nodes[inode,icell]
			coords = d.grid.components[Coordinates]
			node_x = coords[1,inode_k]

			if dims==1
				node_y = 0.0
				n_k = cell_nodes[1,icell]
				n_l = cell_nodes[2,icell]
				x_k = coords[1,n_k]
				x_l = coords[1,n_l]
				node_factor = 0.5*abs(x_l-x_k)
			else
				node_y = coords[2,inode_k]
				node_factor = sys.cellnodefactors[inode,icell]
			end

			idx_vk = idx[1,inode_k]
			idx_nk = idx[2,inode_k]
			idx_pk = idx[3,inode_k]

			v_res = z[idx_nk]-z[idx_pk]
			#v_res -= m.doping(node_x,node_y,region,boundary=false)
			v_res -= m.doping[inode_k]
			v_res *= node_factor*m.q

			ni = m.ni[region]
			tn = m.tn[region]
			tp = m.tp[region]
			recomb = m.recomb(z[idx_nk],z[idx_pk],ni,tn,tp)
			carrier_res = node_factor*recomb

			F[idx_vk] += v_res
			F[idx_nk] += carrier_res
			F[idx_pk] += carrier_res

			j_vk_nk = node_factor*m.q

			J[idx_vk,idx_nk] += j_vk_nk
			J[idx_vk,idx_pk] -= j_vk_nk

			recomb_n(n) = m.recomb(n,z[idx_pk],ni,tn,tp)
			recomb_p(p) = m.recomb(z[idx_nk],p,ni,tn,tp)
			dn_recomb = ForwardDiff.derivative(recomb_n,z[idx_nk])
			dp_recomb = ForwardDiff.derivative(recomb_p,z[idx_pk])

			j_nk_nk = node_factor*dn_recomb
			j_pk_pk = node_factor*dp_recomb

			J[idx_nk,idx_nk] += j_nk_nk
			J[idx_nk,idx_pk] += j_pk_pk
			J[idx_pk,idx_nk] += j_nk_nk
			J[idx_pk,idx_pk] += j_pk_pk

		# for inode in 1: num_nodes(d.grid[CellGeometries][1])
		end

	# for icell in 1: num_cells(d.grid)
	end

	return (nothing)

# function assemble_reaction!
end

# Boundary condition part of residual
function assemble_bcs!(
	F::Zygote.Buffer,
	J::Zygote.Buffer,
	state::NewtonState,
	params::NewtonParams,
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem
)

	m = d.model
	idx = LinearIndices(unknowns(sys))
	dims = dim_space(d.grid)
	z = state.z

	for ibface in 1: num_bfaces(d.grid)

		bregion = d.grid.components[BFaceRegions][ibface]
		bgeom = d.grid.components[BFaceGeometries][1]

		for ibnode in 1: num_nodes(bgeom)

			ibnode_k = d.grid.components[BFaceNodes][ibnode,ibface]

			x = d.grid.components[Coordinates][1,ibnode_k]
			if dims==2
				y = d.grid.components[Coordinates][2,ibnode_k]
			else
				y = 0.0
			end

			idx_vk = idx[1,ibnode_k]
			idx_nk = idx[2,ibnode_k]
			idx_pk = idx[3,ibnode_k]

			if d.b_types[bregion]=="contact"

				ni = m.ni_boundary[bregion]
				#doping = m.doping(x,y,bregion,boundary=true)
				doping = m.doping[ibnode_k]
				n0,p0 = m.doping_to_np(doping,ni)
				v = m.np_to_v(n0,p0,ni,m.vt)

				F[idx_vk] = z[idx_vk]-(v+state.bias[bregion])
				F[idx_nk] = z[idx_nk]-n0
				F[idx_pk] = z[idx_pk]-p0

				F[idx_vk] *= params.dirichlet_scale
				F[idx_nk] *= params.dirichlet_scale
				F[idx_pk] *= params.dirichlet_scale

				for i in 1: num_dof(sys)
					J[idx_vk,i] = 0.0
					J[idx_nk,i] = 0.0
					J[idx_pk,i] = 0.0
				end
				J[idx_vk,idx_vk] = params.dirichlet_scale
				J[idx_nk,idx_nk] = params.dirichlet_scale
				J[idx_pk,idx_pk] = params.dirichlet_scale

			end

		# for ibnode in 1: num_nodes(d.grid[BFaceGeometries][1])
		end

	# for ibface in 1: num_bfaces(d.grid)
	end

	return (nothing)

# function assemble_bcs!
end

# Assemble the residual and Jacobian
function assemble_res_jac(
	state::NewtonState,
	params::NewtonParams,
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem
)

	# The buffer needs to be recreated since it is frozen after each
	# iteration
	# Need to pass freeze=false so that buffer is initialized as zeros
	F = Zygote.Buffer(zero(state.z),false)
	J = Zygote.Buffer(zeros(num_dof(sys),num_dof(sys)),false)

	assemble_flux!(F,J,state,d,sys)
	assemble_reaction!(F,J,state,d,sys)
	assemble_bcs!(F,J,state,params,d,sys)

	# Freeze the buffers and return them as Arrays
	F_vec = copy(F)
	J_mat = sparse(copy(J))
	return (F_vec,J_mat)

# function assemble_res_jac
end

# Assemble flux part of integral
function assemble_integral_flux!(
	state::NewtonState,
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem,
	tf::AbstractVector
)

	m = d.model
	idx = LinearIndices(unknowns(sys))
	dims = dim_space(d.grid)
	z = state.z

	for icell in 1: num_cells(d.grid)

		region = d.grid.components[CellRegions][icell]
		eps = m.e0*m.er[region]
		un = m.mobility_n(0,0,0.0,region,false)
		up = m.mobility_p(0,0,0.0,region,false)

		for iedge in 1: num_edges(d.grid.components[CellGeometries][1])

			if dims==1
				cell_nodes = d.grid.components[CellNodes]
				coords = d.grid.components[Coordinates]
				inode_k = cell_nodes[1,icell]
				inode_l = cell_nodes[2,icell]
				x_k = coords[1,inode_k]
				x_l = coords[1,inode_l]
				edge_factor = 1/abs(x_l-x_k)
			elseif dims==2
				cell_faces = d.grid.components[CellFaces]
				cell_edge = cell_faces[iedge,icell]
				face_nodes = d.grid.components[FaceNodes]
				inode_k = face_nodes[1,cell_edge]
				inode_l = face_nodes[2,cell_edge]
				edge_factor = sys.celledgefactors[iedge,icell]
			end

			idx_vk = idx[1,inode_k]
			idx_vl = idx[1,inode_l]
			idx_nk = idx[2,inode_k]
			idx_nl = idx[2,inode_l]
			idx_pk = idx[3,inode_k]
			idx_pl = idx[3,inode_l]

			deltav = (z[idx_vl]-z[idx_vk])/m.vt
			bp,bm = fbernoulli_pm(deltav)

			n_flux = un*m.vt*(z[idx_nk]*bm-z[idx_nl]*bp)
			p_flux = up*m.vt*(z[idx_pk]*bp-z[idx_pl]*bm)

			edge_current = n_flux-p_flux
			edge_current *= tf[inode_k]-tf[inode_l]
			edge_current *= edge_factor
			state.i += edge_current

		# for iedge in 1: num_edges(d.grid[CellGeometries][1])
		end

	# for icell in 1: num_cells(d.grid)
	end

	return (nothing)

# function assemble_integral_flux!
end

# Integrate to find current
function flux_integrate!(
	state::NewtonState,
	params::NewtonParams,
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem,
	tf::AbstractVector
)

	state.i = 0.0
	assemble_integral_flux!(state,d,sys,tf)
	state.i *= d.model.q
	return (nothing)

# function flux_integrate!
end

# Check if Newton is converged
function check_convergence_flux!(
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
	end
	if state.du/state.du_first < params.tol_rel
		state.converged = status_converged
	end

	# Monotonicity error
	if state.du/state.du_prev > 1/params.tol_mono
		state.converged = err_monotonicity
		error("Monotonicity error")
	end

	# Max iterations reached
	if state.iteration > params.max_iters
		state.converged = err_max_iters
		error("Max iterations reached")
	end

	Zygote.ignore() do
		print_newton(state,params)
	end

	state.iteration += 1
	state.du_prev = state.du
	return (nothing)

# function check_convergence_flux!
end

# Print Newton debug trace
function print_newton(
	state::NewtonState,
	params::NewtonParams
)

	if params.verbose < verbose_print_iters
		return (nothing)
	end

	if state.iteration==1
		printstyled("\nIter\tAbs update\t",color=:light_yellow)
		printstyled("Rel update\tDelta update\t",color=:light_yellow)
		printstyled("Damping\t\tBias\n",color=:light_yellow)
	end

	if state.converged==status_converged
		print_color = :green
	else
		print_color = :default
	end

	status = ""
	status *= @sprintf("%d\t",state.iteration)
	status *= @sprintf("%.5e\t",state.du)
	status *= @sprintf("%.5e\t",state.du/state.du_first)
	status *= @sprintf("%.5e\t",state.du/state.du_prev)
	status *= @sprintf("%.5e\t",state.damp)
	printstyled(status,color=print_color)
	printstyled(state.bias,"\n",color=print_color)

	return (nothing)

# function print_newton
end

# Update damping value
function update_damp!(
	state::NewtonState,
	params::NewtonParams,
	d::Semiconductors.Device
)

	state.damp *= params.damp_growth
	state.damp = min(1.0,state.damp)
	return (nothing)

# function update_damp!
end

# Single Newton step
function newton_step!(
	F::AbstractVector,
	J::AbstractArray,
	state::NewtonState,
	params::NewtonParams,
	d::Semiconductors.Device
)

	state.update_z = sparse(J)\F
	state.du = norm(state.update_z,Inf)

	update_damp!(state,params,d)
	state.z -= state.damp*state.update_z
	return (nothing)

# function newton_step!
end

# Newton solver
function newton!(
	state::NewtonState,
	params::NewtonParams,
	d::Semiconductors.Device,
	sys::VoronoiFVM.AbstractSystem
)

	initialize!(state,params,d,sys)

	while state.converged==status_not_converged
		F,J = assemble_res_jac(state,params,d,sys)
		newton_step!(F,J,state,params,d)
		check_convergence_flux!(state,params)
	end

	return (nothing)

# function newton!
end

