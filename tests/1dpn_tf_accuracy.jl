# 1-D testfunction accuracy test
# Sam Chinnery
# 2022-02-22

module test_1dpn_tf

using VoronoiFVM
using PyPlot
using Zygote
using ExtendableGrids

include("../Semiconductors.jl")
using .Semiconductors

function main()

	d = diode1d(2.0,2.0)

	v0 = equilib(d,damp_initial=0.1)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_equilib = [v0;n0;p0]
	_,sys,ic_bias = non_equilib(d,[0.0 0.0],ic_equilib)
	ic_vec = VoronoiFVM.values(ic_bias)

	factory = VoronoiFVM.TestFunctionFactory(sys)
	tf1 = testfunction(factory,[1],[2])
	tf2 = [0.5+0.25*x for x in vec(d.grid.components[Coordinates])]

	state = NewtonState(ic_vec,[0.0,0.0])
	params = NewtonParams(0.0,0,damp_initial=1.0)

	function solve_current(vf,tf)
		state.bias = [0.0,-vf]
		newton!(state,params,d,sys)
		flux_integrate!(state,params,d,sys,tf)
		return (state.i)
	end

	error_list = []
	vf_range = 0.0: 0.01: 1.0
	for vf in vf_range
		print(vf," ")
		i1 = solve_current(vf,tf1)
		i2 = solve_current(vf,tf2)
		error = abs(i1/i2-1.0)
		push!(error_list,error)
	end
	println()

	semilogy(vf_range,error_list,c="k")
	grid(which="both",axis="y",alpha=0.25)
	grid(alpha=0.25)
	xlabel("Forward voltage (V)")
	ylabel("Relative error")
	tight_layout()

	return (nothing)

end

end

