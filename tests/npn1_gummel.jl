# NPN1 test
# Sam Chinnery
# 2022-02-17

module test_npn1_gummel

using PyPlot
using GridVisualize
using ExtendableGrids

include("../Semiconductors.jl")
using .Semiconductors

function gummel(fldmob)

	default_plotter!(PyPlot)
	d = npn1(5,1,10,5,3,3,0.2,nb=1e5,doping="expbase")
	d.model.fldmob = fldmob

	v0 = equilib(d,damp_initial=0.1,Plotter=nothing)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = [v0;n0;p0]

	vce1 = 0: 0.1: 2.9
	vbe1 = 0.0*ones(length(vce1))
	vee1 = zeros(length(vce1))
	bias_1 = [vee1 vbe1 vce1]
	sol_1 = non_equilib(d,bias_1,ic_bias)[3]

	vbe2 = 0: 0.01: 1.0
	vce2 = 3.0*ones(length(vbe2))
	vee2 = zeros(length(vbe2))
	bias_2 = [vee2 vbe2 vce2]
	j,sys,sol_2 = non_equilib(d,bias_2,sol_1,int_contacts=[2,3],
		damp_initial=0.1)

	return (vbe2,j)

end

function main()

	vbe,j_con = gummel(false)
	vbe,j_fld = gummel(true)

	semilogy(vbe,j_con[: ,1],label="IB, constant mobility",
		c="tab:red",linestyle="-")
	semilogy(vbe,j_con[: ,2],label="IC, constant mobility",
		c="tab:blue",linestyle="-")
	semilogy(vbe,j_fld[: ,1],label="IB, field-dependent mobility",
		c="tab:red",linestyle="--")
	semilogy(vbe,j_fld[: ,2],label="IC, field-dependent mobility",
		c="tab:blue",linestyle="--")

	grid(alpha=0.25)
	grid(which="both",axis="y",alpha=0.25)
	legend()
	xlabel("VBE (V)")
	ylabel("Current (A)")
	title("Gummel plot")
	tight_layout()

	return (nothing)

end

end

