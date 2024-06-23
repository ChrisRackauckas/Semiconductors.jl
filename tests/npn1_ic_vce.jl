# NPN1 I-V using current BC
# Sam Chinnery
# 2022-02-17

module test_npn1_iv

using PyPlot
using VoronoiFVM
using Revise

include("../Semiconductors.jl")
using .Semiconductors

function iv(ib1,fldmob,impact)

	d = npn1(5,1,10,5,3,3,0.2,nb=1e5,doping="expbase",plot_doping=true)
	d.model.fldmob = fldmob
	d.model.impact = impact

	v0 = equilib(d,damp_initial=0.1,Plotter=nothing)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_equilib = [v0;n0;p0]

	vce_list = -1.: 0.1: 5.
	j_list = zeros(length(vce_list))
	vbe = 0.6
	bias = [0.,vbe,vce_list[1]]

	# Setup for BK
	vce_noneq = range(0.,vce_list[1],length=11)
	bias_noneq = hcat(zeros(length(vce_noneq)),vbe*ones(length(vce_noneq)),
		vce_noneq)
	_,sys,ic_bias = non_equilib(d,bias_noneq,ic_equilib)
	ic_bias_vec = VoronoiFVM.values(ic_bias)

	sys.history = VoronoiFVM.DiffEqHistory()
	factory = VoronoiFVM.TestFunctionFactory(sys)
	tf_ib = testfunction(factory,[1,3],[2])
	tf_ic = testfunction(factory,[1,2],[3])

	state = NewtonState(ic_bias_vec,bias)
	params = NewtonParams(ib1,2,damp_growth=1.5,damp_initial=0.1,
		verbose=Semiconductors.verbose_print_iters,max_iters=50,
		do_damp_search=true)

	for i in 1: length(vce_list)
		state.bias[3] = vce_list[i]
		newton_current!(state,params,d,sys,tf_ib)
		printstyled(bias,color=:green)
		println()
		println()
		j_list[i] = current(state.z,d,sys,tf_ic)
	end

	return (vce_list,j_list)

end

function main()

	for ib in [1e-6]
	#for ib in 1e-6: 1e-6: 1e-5
		vce,j1 = iv(ib,false,false)
		plot(vce,j1,label=ib)
	end

	legend()
	xlabel("VCE (V)")
	ylabel("IC (A)")
	title("I-V characteristic")
	tight_layout()

	return (nothing)

end

end

