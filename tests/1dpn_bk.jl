# 1-D diode test using BifurcationKit as solver
# Sam Chinnery
# 2022-01-29

module test_1dpn_bk

using PyPlot
using GridVisualize
using ExtendableGrids
using VoronoiFVM
using BifurcationKit
using LinearAlgebra

include("../Semiconductors.jl")
using .Semiconductors

function iv(fldmob,impact)
	return (nothing)
end

function main()

	d = diode1d(2.0,2.0)

	v0 = equilib(d)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = [v0;n0;p0]

	sys = non_equilib_sys(d,2)
	VoronoiFVM._complete!(sys,create_newtonvectors=true)
	sys.history = VoronoiFVM.DiffEqHistory()

	F(u,p) = bk_residual(sys,u,p)
	J(u,p) = bk_jacobian(sys,u,p)

	optnewton = NewtonPar(tol=1e-10,verbose=true)
	optcont = ContinuationPar(
		dsmin=0.01,
		dsmax=0.2,
		ds=0.1,
		pMin=0.0,
		pMax = 4.1,
		newtonOptions=optnewton
	)
	sol,hist,flag,_ = newton(F,J,ic_bias,[0.0,0.0],optnewton)

	return (norm(sol-ic_bias)/norm(ic_bias))

end

end

