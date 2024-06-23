# VoronoiFVM interface to BifurcationKit
# Sam Chinnery
# 2022-01-29

# Copied from vfvm_diffeq_interface.jl and modified to include parameter
function _bk_res_jac!(sys,u,p)
	uhash = hash((u,p))
	if uhash!=sys.uhash
		ur = reshape(u,sys)
		VoronoiFVM.eval_and_assemble(sys,ur,ur,sys.residual,
			Inf,Inf,0.0,p)
		sys.uhash = uhash
		sys.history.nd += 1
	end
end

# Compute the residual for BifurcationKit.jl (must not be inplace)
function bk_residual(sys,u,p)
	_bk_res_jac!(sys,u,p)
	sys.history.nf += 1
	f = copy(vec(sys.residual))
	return (f)
end

# Compute the Jacobian for BifurcationKit.jl (must not be inplace)
function bk_jacobian(sys,u,p)
	_bk_res_jac!(sys,u,p)
	sys.history.njac += 1
	j = sparse(sys.matrix)
	return (j)
end

