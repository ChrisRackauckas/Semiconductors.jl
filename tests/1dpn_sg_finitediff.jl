# S-G vs finite diff comparison
# Sam Chinnery
# 2022-03-22

module test_1dpn_sg

using PyPlot
using GridVisualize
using ExtendableGrids

include("../Semiconductors.jl")
using .Semiconductors

function main(vf_max)

	d = diode1d(2.0,2.0,ha=0.05,hb=0.01,hc=0.05,na=1e3,nd=1e3)
	x = d.grid.components[Coordinates]

	v0 = equilib(d,damp_initial=0.1)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = [v0;n0;p0]

	vf = 0.0: 0.1: vf_max
	bias_list = [zeros(length(vf)) -vf]
	j,sys,sol = non_equilib(d,bias_list,ic_bias,
		int_contacts=[1],verbose=false)

	# This needs to be run twice, once for each flux discretization
	# The outputs should be saved in the REPL as x, sol_sg and sol_midpt
	return (x,sol)

end

end

