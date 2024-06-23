# 1-D diode test for impact ionization
# Sam Chinnery
# 2022-01-28

module test_1dpn_impact_tf

using PyPlot
using GridVisualize
using ExtendableGrids

include("../Semiconductors.jl")
using .Semiconductors

function iv(fldmob,impact,tf_conc)

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
		verbose=false,max_round=10,catch_conv=true,tf_conc=tf_conc)

	return (vr,j)

end

function main()

	vf,j1 = iv(false,true,false)
	semilogy(vf,j1,label="Standard test function")

	vf,j2 = iv(false,true,true)
	semilogy(vf,j2,label="Concentration-weighted test function")

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

