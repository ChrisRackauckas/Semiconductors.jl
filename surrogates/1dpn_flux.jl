# 1-D diode test with DiffEqFlux interface
# Sam Chinnery
# 2022-03-01

module test_1dpn_flux

using VoronoiFVM
using PyPlot
using Zygote
using ExtendableGrids
using LinearAlgebra
using DiffEqFlux
using GalacticOptim
using BlackBoxOptim

include("../Semiconductors.jl")
using .Semiconductors

function diode_setup(d)

	v0 = equilib(d,damp_initial=0.1)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_equilib = [v0;n0;p0]
	_,sys,ic_bias = non_equilib(d,[0.0 0.0],ic_equilib)
	ic_vec = VoronoiFVM.values(ic_bias)

	tf = [0.5+0.25*x for x in vec(d.grid.components[Coordinates])]

	state = NewtonState(ic_vec,[0.0,0.0])
	params = NewtonParams(0.0,0,
		damp_initial=1.0
	)

	doping_initial = zero(d.grid.components[Coordinates])
	for (i,x) in enumerate(d.grid.components[Coordinates])
		if x<-1e-9
			doping_initial[i] = -1e3
		elseif x>1e-9
			doping_initial[i] = 1e3
		end
	end

	d.model.doping = doping_initial
	return (state,params,sys,tf)

end

function solve_current(state,params,d,sys,tf)
	newton!(state,params,d,sys)
	flux_integrate!(state,params,d,sys,tf)
	return (state.i)
end

function main()

	d_fine = diode1d(2.0,2.0,ha=0.05,hb=0.01,hc=0.05)
	d_coarse = diode1d(2.0,2.0,ha=0.8,hb=0.16,hc=0.8)
	d_coarse_init = diode1d(2.0,2.0,ha=0.8,hb=0.16,hc=0.8)

	state_fine,params_fine,sys_fine,tf_fine = diode_setup(d_fine)
	state_coarse,params_coarse,sys_coarse,tf_coarse = diode_setup(d_coarse)
	state_init,params_init,sys_init,tf_init = diode_setup(d_coarse_init)
	params_coarse.max_iters = 1000

	coarse_doping_initial = vec(d_coarse.model.doping)

	i_list_fine = []
	vf_list = 0.1: 0.1: 1.0
	for vf in vf_list
		print(vf," ")
		state_fine.bias = [0.0,-vf]
		i_fine = solve_current(
			state_fine,
			params_fine,
			d_fine,
			sys_fine,
			tf_fine
		)
		push!(i_list_fine,i_fine)
	end
	println()

	function loss(theta)

		weight = theta[1]
		weight_sig = 1.0/(1.0+exp(-weight))

		#d_coarse.model.doping = coarse_doping_initial
		#doping_residual = 1.0 .+ weight_sig.*tanh.(theta[2: end])
		#doping = vec(d_coarse.model.doping).*doping_residual
		d_coarse.model.doping = theta[2: end]

		charge_neutrality_ic!(state_coarse,d_coarse,sys_coarse)
		error_total = 0.0

		for i in 1: length(vf_list)
			vf = vf_list[i]
			state_coarse.bias = [0.0,-vf]
			i_coarse = solve_current(
				state_coarse,
				params_coarse,
				d_coarse,
				sys_coarse,
				tf_coarse
			)
			error = i_coarse/i_list_fine[i]-1.0
			error_total += error*error
		end

		error_total = sqrt(error_total)
		return (error_total)

	end

	loss_accum = []
	theta_accum = []
	iteration = 1
	max_iters = 10

	function cb(theta,loss)
		append!(loss_accum,loss)
		push!(theta_accum,theta)
		println("Iteration ",iteration)
		println("\tLoss: ",loss)
		#println("\ttheta: ",theta)
		println()
		iteration += 1
		return (false)
	end

	theta_initial = zeros(1+length(d_coarse.grid[Coordinates]))
	theta_initial[2: end] .= coarse_doping_initial
	#theta_initial[1] = -log(90.0)
	loss_initial = loss(theta_initial)
	println()
	println("Initial loss: ",loss_initial)
	println()

	high_fact = 2.0
	low_fact = 1/high_fact
	lb = [t>=0 ? t*low_fact : t*high_fact for t in theta_initial]
	ub = [t>=0 ? t*high_fact : t*low_fact for t in theta_initial]

	res = DiffEqFlux.sciml_train(
		loss,
		theta_initial,
		BFGS(),
		#ADAM(0.1),
		#BBO_adaptive_de_rand_1_bin_radiuslimited(),
		GalacticOptim.AutoZygote(),
		#lower_bounds=lb,
		#upper_bounds=ub,
		cb=cb,
		maxiters=max_iters
	)

	d_fine = diode1d(2.0,2.0,ha=0.05,hb=0.01,hc=0.05)
	d_coarse_init = diode1d(2.0,2.0,ha=0.8,hb=0.16,hc=0.8)

	theta_min_idx = findmin(loss_accum)[2]
	theta = theta_accum[theta_min_idx]
	weight = theta[1]
	weight_sig = 1.0/(1.0+exp(-weight))

	#d_coarse.model.doping = vec(coarse_doping_initial)
	#doping_residual = 1.0 .+ weight_sig.*tanh.(theta[2: end])
	#d_coarse.model.doping = vec(d_coarse.model.doping).*doping_residual
	d_coarse.model.doping .= theta[2: end]

	state_fine,params_fine,sys_fine,tf_fine = diode_setup(d_fine)
	state_init,params_init,sys_init,tf_init = diode_setup(d_coarse_init)

	charge_neutrality_ic!(state_coarse,d_coarse,sys_coarse)

	i_list_fine = []
	vf_list_test = 0.01: 0.01: 1.00
	error_initial = []
	error_optimized = []
	for i in 1: length(vf_list_test)
		vf = vf_list_test[i]
		state_coarse.bias = [0.0,-vf]
		state_init.bias = [0.0,-vf]
		state_fine.bias = [0.0,-vf]

		i_fine = solve_current(
			state_fine,
			params_fine,
			d_fine,
			sys_fine,
			tf_fine
		)
		push!(i_list_fine,i_fine)

		i_coarse = solve_current(
			state_coarse,
			params_coarse,
			d_coarse,
			sys_coarse,
			tf_coarse
		)
		error = i_coarse/i_list_fine[i]-1.0
		push!(error_optimized,error)

		i_initial = solve_current(
			state_init,
			params_init,
			d_coarse_init,
			sys_init,
			tf_init
		)
		error = i_initial/i_list_fine[i]-1.0
		push!(error_initial,error)
	end

	plot(vf_list_test,error_initial,c="k",linestyle="--",label="Initial")
	plot(vf_list_test,error_optimized,c="k",label="Optimized")
	xlabel("Forward voltage (V)")
	ylabel("Relative error")
	legend()
	tight_layout()

	return (loss_accum,theta_accum,coarse_doping_initial,theta[2: end])

end

end

