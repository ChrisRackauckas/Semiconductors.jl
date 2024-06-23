# 1-D diode test for impact ionization
# Sam Chinnery
# 2022-01-28

module test_1dpn_impact

using PyPlot
using GridVisualize
using ExtendableGrids

include("../Semiconductors.jl")
using .Semiconductors

function iv(fldmob,impact)

	default_plotter!(PyPlot)
	d = diode1d(2.0,2.0,ha=1e-5,hb=0.1,hc=1e-5)
	d.model.fldmob = fldmob
	d.model.impact = impact

	v0 = equilib(d)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = [v0;n0;p0]

	vr = 0.0: 0.1: 120.0
	bias_list = [zeros(length(vr)) vr]
	j,sys,sol = non_equilib(d,bias_list,ic_bias,int_contacts=[2],
		verbose=false,max_round=10,catch_conv=true)

	return (vr,j)

end

function main()

	vf,j1 = iv(false,false)
	semilogy(vf,j1,label="No field, no impact")

	vf,j2 = iv(true,false)
	semilogy(vf,j2,label="Field, no impact")

	vf,j3 = iv(false,true)
	semilogy(vf,j3,label="No field, impact")

	vf,j4 = iv(true,true)
	semilogy(vf,j4,label="Field, impact")

	grid(alpha=0.25)
	grid(which="both",axis="y",alpha=0.25)
	legend()
	xlabel("Reverse bias voltage (V)")
	ylabel("Current (A/um2)")
	title("I-V characteristics, 1-dimensional diode")
	tight_layout()
	return (nothing)

end

end

