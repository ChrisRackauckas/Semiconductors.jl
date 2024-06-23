# Current BC test
# Sam Chinnery
# 2022-02-10

module test_current_bc

using VoronoiFVM
using PyPlot

include("../Semiconductors.jl")
using .Semiconductors

function main()

	d = diode1d(2.0,2.0)
	d.model.fldmob = false
	d.model.impact = true

	v0 = equilib(d)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = [v0;n0;p0]

	vr = vcat(1.: 2.: 113.)
	#vr = [0.]
	bias_list = hcat(zeros(length(vr),1),vr)
	j,sys,sol = non_equilib(d,bias_list,ic_bias,int_contacts=2)
	semilogy(vr,j,c="k",marker="o",mfc="k",mec="None")

	sys.history = VoronoiFVM.DiffEqHistory()
	factory = VoronoiFVM.TestFunctionFactory(sys)
	tf = testfunction(factory,[1],[2])

	j_initial = j[end]
	j_step = 1e-10
	j_length = 11
	j_range = range(j_initial+j_step,step=j_step,length=j_length)
	println(j_range)
	println()

	z = VoronoiFVM.values(sol)
	bias = bias_list[end,: ]

	for j_set in j_range
		z,bias = newton_current(j_set,d,sys,tf,z,bias,2)
	end

	return (nothing)

end

end

