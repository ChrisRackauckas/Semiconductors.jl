# 1-D diode test using BifurcationKit to plot continuation
# Sam Chinnery
# 2022-02-02

module test_1dpn_bk_iv

using GridVisualize
using ExtendableGrids
using VoronoiFVM
using BifurcationKit
using Plots
using Setfield
using LinearAlgebra

include("../Semiconductors.jl")
using .Semiconductors

function main()

	VoronoiFVM.set_dirichlet(1e30)

	d = diode1d(2.0,2.0,ha=1e-5,hb=0.1,hc=1e-5)
	d.model.fldmob = false
	d.model.impact = false

	v0 = equilib(d)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = VoronoiFVM.values([v0;n0;p0])

	sys = non_equilib_sys(d,2)
	VoronoiFVM._complete!(sys,create_newtonvectors=true)
	sys.history = VoronoiFVM.DiffEqHistory()
	idx = unknown_indices(unknowns(sys))

	factory = VoronoiFVM.TestFunctionFactory(sys)
	tf = testfunction(factory,[2],[1])
	function current(u,p)
		ur = reshape(u,sys)
		bflux = do_integrate(sys,tf,ur,d,idx)
		j = d.model.q*(bflux[3]-bflux[2])
		return (log10(abs(j)))
	end

	F(u,p) = bk_residual(sys,u,p)
	J(u,p) = bk_jacobian(sys,u,p)

	optnewton = NewtonPar(tol=1e-10,verbose=true)
	sol,hist,flag,_ = newton(F,J,ic_bias,[0.0,1.0],optnewton)

	function plot_scalar(u,p;k...)
		ur = reshape(u,sys)
		plot!(vec(d.grid[Coordinates]),log10.(vec(ur[2,: ]));
			xlabel="Position (um)",ylabel="Carrier concentration",
			label="",k...)
		plot!(vec(d.grid[Coordinates]),log10.(vec(ur[3,: ]));
			xlabel="Position (um)",ylabel="Carrier concentration",
			label="",k...)
		return (nothing)
	end

	function dot_current(u1,u2)
		u1r = reshape(u1,sys)
		u2r = reshape(u2,sys)
		bflux1 = do_integrate(sys,tf,u1r,d,idx)
		bflux2 = do_integrate(sys,tf,u2r,d,idx)
		j1 = d.model.q*(bflux1[3]-bflux1[2])
		j2 = d.model.q*(bflux2[3]-bflux2[2])
		dot_j = abs(log10(max(0,j1))*log10(max(0,j2)))
		dot_j = dot_j==Inf ? 0 : dot_j
		#println(j1,", ",j2,", ",dot_j)
		return (dot_j)
	end
	dt = BifurcationKit.DotTheta(dot_current)

	optnewton2 = NewtonPar(
		tol=1e-10,
		verbose=true,
		maxIter=25,
		linesearch=true,
		eigsolver=EigArpack()
	)
	optcont = ContinuationPar(
		ds=10.0,
		dsmin=1.0,
		dsmax=1000.0,
		pMin=0.1,
		pMax=200.0,
		theta=0.5,
		doArcLengthScaling=false,
		maxSteps=50,
		a=0.5,
		detectBifurcation=0,
		newtonOptions=optnewton2
	)
	br, = continuation(F,J,sol,[0.0,1.0],(@lens _[2]),optcont;
		plot=true,
		plotSolution=plot_scalar,
		verbosity=2,
		recordFromSolution=current,
		tangentAlgo=SecantPred()
		#tangentAlgo=PolynomialPred(SecantPred(),2,5,ic_bias)
		#tangentAlgo=MultiplePred(copy(ic_bias),0.01,20),
		#dotPALC=dt
	)

	return (nothing)

end

end

