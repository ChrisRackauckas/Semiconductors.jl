# Deflated test
# Sam Chinnery
# 2022-02-09

module test_deflation

using VoronoiFVM
using BifurcationKit
using PyPlot
using LinearAlgebra

include("../Semiconductors.jl")
using .Semiconductors

function main(p,a)

	VoronoiFVM.set_dirichlet(1e30)

	d = diode1d(2.0,2.0,ha=1e-5,hb=0.1,hc=1e-5)
	d.model.fldmob = false
	d.model.impact = true

	v0 = equilib(d)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = [v0;n0;p0]

	vr_init = 100.
	vr = 1.: 4.: vr_init
	bias_list = [zeros(length(vr)) vr]
	j,sys,ic_noneq = non_equilib(d,bias_list,ic_bias,int_contacts=2)
	plot(vr,j,c="k")
	ic_noneq_vec = VoronoiFVM.values(ic_noneq)

	sys.history = VoronoiFVM.DiffEqHistory()
	factory = VoronoiFVM.TestFunctionFactory(sys)
	tf = testfunction(factory,[2],[1])
	idx = unknown_indices(unknowns(sys))

	optnewton = NewtonPar(
		verbose=true,
		tol=1e-10,
		maxIter=200
	)
	F(u,p) = bk_residual(sys,u,p)
	J(u,p) = bk_jacobian(sys,u,p)

	v_bk_initial = 105.
	p_bk_initial = [0.,v_bk_initial]
	sol_bk_initial, = newton(F,J,ic_noneq_vec,p_bk_initial,optnewton)

	sol_bk_initial_r = reshape(sol_bk_initial,sys)
	bflux1 = do_integrate(sys,tf,sol_bk_initial_r,d,idx)
	j1 = d.model.q*(bflux1[3]-bflux1[2])
	println(j1)

	norminf(x) = norm(x,Inf)
	optnewton_defl = NewtonPar(
		verbose=false,
		tol=1e-10,
		maxIter=200
	)
	deflation_op = DeflationOperator(p,dot,a,[sol_bk_initial])

	vr_defl = 105: -1.: 80.
	bias_list_defl = [zeros(length(vr_defl)) vr_defl]
	for i in 1: 100
		sol_bk_defl, = newton(
			F,J,
			#vec(ic_noneq_vec.*(0.995.+rand(num_dof(sys),1)*0.01)),
			ic_noneq_vec,
			p_bk_initial,
			optnewton_defl,
			deflation_op,
			normN=norminf
		)
		sol_bk_defl_r = reshape(sol_bk_defl,sys)
		bflux2 = do_integrate(sys,tf,sol_bk_defl_r,d,idx)
		j2 = d.model.q*(bflux2[3]-bflux2[2])
		println()
		println(j2)
		push!(deflation_op,sol_bk_defl)
		j_list, = non_equilib(
			d,bias_list_defl,sol_bk_defl_r,int_contacts=2)
		#println(j_list)
		if j2<0
			plot(vr_defl,j_list,c="tab:red")
		else
			plot(vr_defl,j_list,c="tab:blue")
		end
	end

	return (nothing)

end

end

