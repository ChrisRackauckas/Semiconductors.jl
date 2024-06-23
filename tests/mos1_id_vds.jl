# MOS1 ID-VGS test
# Sam Chinnery
# 2022-01-28

module test_mos1_id_vds

using PyPlot
using GridVisualize
using ExtendableGrids

include("../Semiconductors.jl")
using .Semiconductors

function id_vds(vgs,fldmob,impact)

	d = mos1(0.1,0.05,0.2,0.05,0.002,0.1,0.025)
	d.model.fldmob = fldmob
	d.model.impact = impact

	v0 = equilib(d,verbose=true)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = [v0;n0;p0]

	vgs1 = 0.0: 0.1: vgs
	vds1 = 0.0*ones(length(vgs1))
	vbs1 = zeros(length(vgs1))
	vss1 = zeros(length(vgs1))
	bias_1 = [vgs1 vss1 vds1 vbs1]
	sol_1 = non_equilib(d,bias_1,ic_bias,verbose=true,max_round=10)[3]

	vds2 = 0.0: 0.01: 2.0
	vgs2 = vgs*ones(length(vds2))
	vbs2 = zeros(length(vds2))
	vss2 = zeros(length(vds2))
	bias_2 = [vgs2 vss2 vds2 vbs2]
	j,sys,sol_2 = non_equilib(d,bias_2,sol_1,int_contacts=[3],
		verbose=true,Plotter=nothing,catch_conv=true,max_round=10)

	if !impact
		plot(vds2,j,label=vgs,c="k")
	else
		plot(vds2,j,label=vgs,c="k",linestyle="--")
	end
	return (nothing)

end

function main()

	default_plotter!(PyPlot)

	for vgs in 0.0: 0.1: 1.0
		id_vds(vgs,true,false)
		id_vds(vgs,true,true)
	end

	grid(alpha=0.25)
	grid(which="both",axis="y",alpha=0.25)
	xlabel("VDS (V)")
	ylabel("ID (A/um)")
	title("I-V characteristic")
	tight_layout()

	return (nothing)

end

end

