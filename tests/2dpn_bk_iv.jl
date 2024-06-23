# 1-D diode test using BifurcationKit to plot continuation
# Sam Chinnery
# 2022-02-01

module test_2dpn_bk_iv

using GridVisualize
using ExtendableGrids
using VoronoiFVM
using BifurcationKit
using Plots
using Setfield

include("../Semiconductors.jl")
using .Semiconductors

function main()

	VoronoiFVM.set_dirichlet(1e30)

	d = diode2d(2.0,2.0,1.0,0.5)
	d.model.fldmob = false
	d.model.impact = true

	v0 = equilib(d)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = VoronoiFVM.values([v0;n0;p0])

	sys = non_equilib_sys(d,2)
	VoronoiFVM._complete!(sys,create_newtonvectors=true)
	sys.history = VoronoiFVM.DiffEqHistory()
	idx = unknown_indices(unknowns(sys))

	# Use VFVM integration method if generic operator is not needed
	function do_integrate(sys,tf,sol,device,idx)
		if d.model.fldmob || d.model.impact return (
			Semiconductors.integrate_generic(sys,tf,sol,device,idx))
		else return (VoronoiFVM.integrate_stdy(sys,tf,sol)) end
	end

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
		plot!(vec(d.grid[Coordinates]),vec(ur[1,: ]);
			xlabel="Position (um)",ylabel="Potential (V)",
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
		dot_j = j1*j2
		#println(dot_j)
		return (dot_j)
	end
	dt = BifurcationKit.DotTheta(dot_current)

	optnewton2 = NewtonPar(tol=1e-10,verbose=true,maxIter=25,
		linesearch=true)
	optcont = ContinuationPar(
		ds=0.1,
		dsmin=1e-4,
		dsmax=10.0,
		pMin=1.0,
		pMax=120.0,
		theta=0.0001,
		doArcLengthScaling=false,
		maxSteps=100,
		a=0.5,
		newtonOptions=optnewton2
	)
	br, = continuation(F,J,sol,[0.0,1.0],(@lens _[2]),optcont;
		plot=true,
		verbosity=3,
		recordFromSolution=current,
		tangentAlgo=SecantPred(),
		#dotPALC=dt
	)

	return (nothing)

end

end

