# 1-D diode test with Newton PALC
# Sam Chinnery
# 2022-02-05

module test_continuation_newton

using VoronoiFVM
using PyPlot

include("../Semiconductors.jl")
using .Semiconductors

function main()

	VoronoiFVM.set_dirichlet(1e30)

	d = diode1d(2.0,2.0)
	d.model.fldmob = false
	d.model.impact = true

	v0 = equilib(d)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = [v0;n0;p0]

	vr_init = 1.
	vr = [1.]
	bias_list = [zeros(length(vr)) vr]
	_,sys,ic_noneq = non_equilib(d,bias_list,ic_bias)
	ic_noneq_vec = VoronoiFVM.values(ic_noneq)

	sys.history = VoronoiFVM.DiffEqHistory()
	factory = VoronoiFVM.TestFunctionFactory(sys)
	tf = testfunction(factory,[2],[1])

	state = IVState(ic_noneq_vec,[0.,vr_init],2,1.)
	curve = IVCurve(state,0.,10.,1e-5)

	state.scale_j = 1e10
	state.max_iters = 25
	#state.theta = 1.99
	#state.log_j = true
	#state.damp_initial = 0.1
	#state.n_damp = 3
	#state.use_n1 = false
	state.verbose = Semiconductors.verbose_print_cont

	curve.max_steps = 2000
	curve.n_solve_max = 5000
	curve.tol_trunc = 1.

	continuation!(state,curve,d,sys,tf)
	semilogy(curve.v,abs.(curve.i),
		linestyle="None",marker="o",mfc="k",mec="None")

	return (nothing)

end

end

