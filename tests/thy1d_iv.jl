# 1-D thyristor test
# Sam Chinnery
# 2022-03-08

module test_thy1d_iv

using PyPlot
using GridVisualize
using ExtendableGrids

include("../Semiconductors.jl")
using .Semiconductors

function iv()

	# Modified to be a PN diode with heavily doped end regions
	d = thy1d(x1=-3.0,x2=-2.0,x3=0.0,x4=2.0,x5=3.0,
		h1=0.2,h2=0.02,h3=0.1,h4=0.02,h5=0.2,
		n1=1.0e6,n2=-1.0e3,n3=-1.0e3,n4=1.0e6)
	d.model.impact = true

	v0 = equilib(d,damp_initial=0.1)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = [v0;n0;p0]

	vf = 0.0: 0.1: 1.0
	bias_list = [zeros(length(vf)) vf]
	j,sys,sol = non_equilib(d,bias_list,ic_bias,
		int_contacts=[2],verbose=false)

	#semilogy(vf,j,c="k")
	#grid(which="both",axis="y",alpha=0.25)
	#grid(alpha=0.25)
	return (nothing)

end

function main()
	iv()
	return (nothing)
end

end

