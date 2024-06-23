# 2-D diode test
# Sam Chinnery
# 2022-01-20

module test_2dpn_iv

using PyPlot
using GridVisualize
using ExtendableGrids

include("../Semiconductors.jl")
using .Semiconductors

function iv(fldmob,bgmesh)

	default_plotter!(PyPlot)
	d = diode2d(2.0,2.0,1.0,0.5,bgmesh=bgmesh)
	d.model.fldmob = fldmob

	v0 = equilib(d,damp_initial=0.1,Plotter=nothing)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = [v0;n0;p0]

	vf = 0: 1.0: 90.0
	bias_list = [zeros(length(vf)) vf]
	j,sys,sol = non_equilib(d,bias_list,ic_bias,
		int_contacts=[1],verbose=true,Plotter=nothing)

	return (vf,j)

end

function main()

	vf,j_con_nobg = iv(false,false)
	vf,j_fld_nobg = iv(true,false)
	vf,j_con_bg = iv(false,true)
	vf,j_fld_bg = iv(true,true)

	semilogy(vf,-j_con_nobg,
		label="Constant mobility, no background mesh")
	semilogy(vf,-j_fld_nobg,
		label="Field-dependent mobility, no background mesh")
	semilogy(vf,-j_con_bg,
		label="Constant mobility, background mesh")
	semilogy(vf,-j_fld_bg,
		label="Field-dependent mobility, background mesh")

	grid(alpha=0.25)
	grid(which="both",axis="y",alpha=0.25)
	ylim(1e-13,1e-5)
	legend()
	xlabel("Forward bias voltage (V)")
	ylabel("Current (A)")
	title("I-V characteristics, 2-dimensional diode")
	tight_layout()

	return (nothing)

end

end

