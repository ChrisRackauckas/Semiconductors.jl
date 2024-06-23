# NPN1 Beta-IC curve
# Sam Chinnery
# 2022-02-17

module test_npn1_beta_ic

using PyPlot
using VoronoiFVM

include("../Semiconductors.jl")
using .Semiconductors

function curve(vce,fldmob,impact)

	d = npn1(5,1,10,5,3,3,0.2,nb=1e5,doping="expbase")
	d.model.fldmob = fldmob
	d.model.impact = impact

	v0 = equilib(d,damp_initial=0.1,Plotter=nothing)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_equilib = [v0;n0;p0]

	vbe = 0.3
	vce_noneq = range(0.,vce,length=11)
	bias_noneq = hcat(zeros(length(vce_noneq)),vbe*ones(length(vce_noneq)),
		vce_noneq)
	_,sys,ic_bias = non_equilib(d,bias_noneq,ic_equilib)
	ic_bias_vec = VoronoiFVM.values(ic_bias)

	sys.history = VoronoiFVM.DiffEqHistory()
	factory = VoronoiFVM.TestFunctionFactory(sys)
	tf_ib = testfunction(factory,[1,3],[2])
	tf_ic = testfunction(factory,[1,2],[3])

	ic_list = 10.0.^range(-12,-4,length=51)
	ib_list = zero(ic_list)
	bias = [0.,vbe,vce]

	state = NewtonState(ic_bias_vec,bias)
	params = NewtonParams(0.0,2,damp_growth=1.5,damp_initial=0.1,
		verbose=Semiconductors.verbose_print_iters,max_iters=50,
		do_damp_search=false)

	for i in 1: length(ic_list)
		params.set_current = ic_list[i]
		newton_current!(state,params,d,sys,tf_ic)
		printstyled(bias,color=:green)
		println()
		println()
		ib_list[i] = current(state.z,d,sys,tf_ib)
	end

	return (ic_list,ib_list)

end

function main()

	ic,ib = curve(3.0,false,false)
	semilogx(ic,ic./ib,c="k",marker="o",mfc="k",mec="None")
	grid(which="both",axis="x",alpha=0.25)
	grid(alpha=0.25)

	xlabel("IC (A)")
	ylabel("Beta")
	title("DC Current gain")
	tight_layout()

	return (nothing)

end

end

