# MOS1 ID-VGS test
# Sam Chinnery
# 2022-01-26

module test_mos1_id_vgs

using PyPlot
using GridVisualize
using ExtendableGrids

include("../Semiconductors.jl")
using .Semiconductors

function id_vgs(fldmob)

	default_plotter!(PyPlot)
	d = mos1(0.1,0.05,0.2,0.05,0.002,0.1,0.025)
	d.model.fldmob = fldmob

	v0 = equilib(d,verbose=true)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = [v0;n0;p0]

	vds1 = 0.1: 0.1: 0.9
	vgs1 = 0.0*ones(length(vds1))
	vbs1 = zeros(length(vds1))
	vss1 = zeros(length(vds1))
	bias_1 = [vgs1 vss1 vds1 vbs1]
	sol_1 = non_equilib(d,bias_1,ic_bias,verbose=true)[3]

	vgs2 = 0.0: 0.01: 1.0
	vds2 = 1.0*ones(length(vgs2))
	vbs2 = zeros(length(vgs2))
	vss2 = zeros(length(vgs2))
	bias_2 = [vgs2 vss2 vds2 vbs2]
	j,sys,sol_2 = non_equilib(d,bias_2,sol_1,int_contacts=[3],
		verbose=true,Plotter=nothing)

	return (vgs2,j)

end

function main()

	vgs,j_con = id_vgs(false)
	vgs,j_fld = id_vgs(true)

	semilogy(vgs,j_con,label="Constant mobility",c="k",linestyle="-")
	semilogy(vgs,j_fld,label="Field-dependent mobility",c="k",
		linestyle="--")

	grid(alpha=0.25)
	grid(which="both",axis="y",alpha=0.25)
	legend()
	xlabel("VGS (V)")
	ylabel("ID (A/um)")
	title("I-V characteristic")
	tight_layout()

	return (nothing)

end

end

