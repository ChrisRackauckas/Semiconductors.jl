# 1-D diode test
# Sam Chinnery
# 2022-01-20

module test_1dpn_iv

using PyPlot
using GridVisualize
using ExtendableGrids

include("../Semiconductors.jl")
using .Semiconductors

function iv(fldmob)

	default_plotter!(PyPlot)
	d = diode1d(2.0,2.0,ha=0.05,hb=0.01,hc=0.05,na=1e3,nd=1e3)
	d.model.fldmob = fldmob

	v0 = equilib(d,damp_initial=0.1,Plotter=nothing)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = [v0;n0;p0]

	vf = 0.0: 0.01: 1.0
	bias_list = [zeros(length(vf)) -vf]
	j,sys,sol = non_equilib(d,bias_list,ic_bias,
		int_contacts=[1],verbose=false)

	if d.model.fldmob
		semilogy(vf,j,c="k",linestyle="--",label="Field-dependent \
		mobility")
	else
		semilogy(vf,j,c="k",linestyle="-",label="Constant mobility")
	end

	return (nothing)

end

function main()
	iv(false)
	iv(true)
	grid(alpha=0.25)
	grid(which="both",axis="y",alpha=0.25)
	legend()
	xlabel("Forward bias voltage (V)")
	ylabel("Current (A)")
	title("I-V characteristics, 1-dimensional diode")
	tight_layout()
	return (nothing)
end

end

