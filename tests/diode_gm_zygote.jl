# Diode transconductance plot via Zygote AD
# Sam Chinnery
# 2022-02-22

module test_diode_gm

using VoronoiFVM
using PyPlot
using Zygote

include("../Semiconductors.jl")
using .Semiconductors

function main()

	#d = diode1d(2.0,2.0)
	d = diode2d(2.0,2.0,1.0,0.5)

	v0 = equilib(d,damp_initial=0.1)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_equilib = [v0;n0;p0]
	_,sys,ic_bias = non_equilib(d,[0.0 0.0],ic_equilib)
	ic_vec = VoronoiFVM.values(ic_bias)

	factory = VoronoiFVM.TestFunctionFactory(sys)
	tf = testfunction(factory,[1],[2])

	state = NewtonState(ic_vec,[0.0,0.0])
	params = NewtonParams(0.0,0,damp_initial=1.0,
		verbose=verbose_print_iters)

	function solve_current(vf)
		state.bias = [0.0,-vf]
		newton!(state,params,d,sys)
		flux_integrate!(state,params,d,sys,tf)
		return (state.i)
	end

	vf_list = 0.0: 0.1: 1.0
	i_list = []
	gm_list = []
	for vf in vf_list
		i,gm = withgradient(solve_current,vf)
		push!(i_list,i)
		push!(gm_list,gm)
	end

	plot(vf_list,i_list,c="k",label="Current")
	plot(vf_list,gm_list,c="k",linestyle="--",label="Transconductance")
	legend()
	tight_layout()

	return (nothing)

end

end

