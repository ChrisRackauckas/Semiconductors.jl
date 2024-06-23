# 1-D diode test function comparison
# Sam Chinnery
# 2022-04-01

module test_1dpn_tf

using PyPlot
using ExtendableGrids
using VoronoiFVM

include("../Semiconductors.jl")
using .Semiconductors

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["text.latex.preamble"] = "\\usepackage{pxfonts}\n\\usepackage{siunitx}"

rcParams["figure.figsize"] = [6.0, 3.0]

function iv()

	d = diode1d(2.0,2.0)

	v0 = equilib(d,damp_initial=0.1,Plotter=nothing)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = [v0;n0;p0]

	vf = 0.0: 0.1: 1.0
	bias_list = [zeros(length(vf)) -vf]
	j,sys,sol,tf_n_fwd,tf_p_fwd = non_equilib(d,bias_list,ic_bias,
		int_contacts=[1],tf_conc=true,return_tfs=true)

	vr = 0.0: 0.1: 1.0
	bias_list = [zeros(length(vf)) vr]
	j,sys,sol,tf_n_rev,tf_p_rev = non_equilib(d,bias_list,ic_bias,
		int_contacts=[1],tf_conc=true,return_tfs=true)

	return (d,j,sol,tf_n_fwd[1],tf_p_fwd[1],tf_n_rev[1],tf_p_rev[1])

end

function main()

	d,j,sol,tf_n_fwd,tf_p_fwd,tf_n_rev,tf_p_rev = iv()
	x = vec(d.grid.components[Coordinates])
	v = vec(sol[1,: ])
	n = vec(sol[2,: ])
	p = vec(sol[3,: ])
	un = d.model.mobility_n(0.0,0.0,0.0,0.0,false)
	up = d.model.mobility_p(0.0,0.0,0.0,0.0,false)
	tf_laplace = 0.25*(2.0.-x)

	i = 0.0
	i1 = 0.0
	for k in 1: length(n)-1
		dv = (v[k+1]-v[k])/d.model.vt
		bp,bm = fbernoulli_pm(dv)
		gn = un*(bm*n[k]-bp*n[k+1])
		gp = up*(bp*p[k]-bm*p[k+1])
		i += gn-gp
		if k==1
			i1 = d.model.q*d.model.vt*(gn-gp)/(x[2]-x[1])
		end
	end
	i *= 0.25*d.model.q*d.model.vt

	fig,axs = subplots(1,2)

	axs[1].plot(x,vec(tf_n_fwd),c="k",label="Electrons")
	axs[1].plot(x,vec(tf_p_fwd),c="k",linestyle="--",label="Holes")
	#axs[1].plot(x,tf_laplace,c="k",linestyle=":",label="Laplace")
	axs[1].legend()
	axs[1].set_xlabel("Position \$(\\SI{}{\\micro m})\$")
	axs[1].set_ylabel("\$T_1(x)\$")
	axs[1].grid(alpha=0.25)
	axs[1].set_title("Forward-biased \$(V_F=\\SI{1}{V})\$")

	axs[2].plot(x,vec(tf_n_rev),c="k",label="Electrons")
	axs[2].plot(x,vec(tf_p_rev),c="k",linestyle="--",label="Holes")
	#axs[2].plot(x,tf_laplace,c="k",linestyle=":",label="Laplace")
	axs[2].legend()
	axs[2].set_xlabel("Position \$(\\SI{}{\\micro m})\$")
	axs[2].set_ylabel("\$T_1(x)\$")
	axs[2].grid(alpha=0.25)
	axs[2].set_title("Reverse-biased \$(V_R=\\SI{1}{V})\$")

	fig.tight_layout()
	savefig("nanz_tf_1d.pdf")

	println("Test function: ",j[end])
	println("Direct: ",i)
	println("First grid point: ",i1)
	return (nothing)

end

end

