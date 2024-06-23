# Deflated test
# Sam Chinnery
# 2022-02-14

module test_deflation

using VoronoiFVM
using BifurcationKit
using PyPlot
using LinearAlgebra

include("../Semiconductors.jl")
using .Semiconductors

function main(p,a)

	# Device setup
	d = diode1d(2.0,2.0)
	d.model.fldmob = false
	d.model.impact = true

	# Equilibrium ICs
	v0 = equilib(d)
	n0 = d.model.n_dist.(v0,0,d.model.ni[1],d.model.vt)
	p0 = d.model.p_dist.(v0,0,d.model.ni[1],d.model.vt)
	ic_bias = [v0;n0;p0]

	# Create non-equilibrium system
	_,sys,ic_main = non_equilib(d,[0. 0.],ic_bias,tol_absolute=1e-13)
	ic_main = VoronoiFVM.values(ic_main)

	# Voltages to trace along main branch
	vr_main = 1.: 1.: 107.
	bias_main = hcat(zeros(length(vr_main)),vr_main)

	# Plot the main branch
	j,_,_ = non_equilib(d,bias_main,ic_bias,int_contacts=2)
	plot(vr_main,j,c="k")

	# Offsets to trace at each deflation branch
	vr_deflation = 0.: -1.: -10.

	# Number of deflation attempts at each voltage on the main branch
	n_deflation = 10

	# Main branch and deflation branch I-V curves
	i_main = zeros(length(vr_main))
	v_main = zeros(length(vr_main))
	i_deflation = zeros(length(vr_deflation))
	v_deflation = zeros(length(vr_deflation))

	# Newton options for deflation
	optnewton_defl = NewtonPar(
		verbose=false,
		tol=1e-13,
		maxIter=200
	)

	# Setup for BK
	sys.history = VoronoiFVM.DiffEqHistory()
	factory = VoronoiFVM.TestFunctionFactory(sys)
	tf = testfunction(factory,[2],[1])
	idx = unknown_indices(unknowns(sys))
	F(u,p) = bk_residual(sys,u,p)
	J(u,p) = bk_jacobian(sys,u,p)

	# Loop through main branch voltages
	for i in 1: length(vr_main)

		# Advance along the main branch
		bias = bias_main[i,: ]
		println(bias)
		ic_main, = newton(F,J,ic_main,bias,optnewton_defl)

		# Look for another branch and continue it
		deflation_op = DeflationOperator(p,dot,a,[ic_main])
		for n in 1: n_deflation
			try
				ic_main_p = vec(
					ic_main.*(0.995.+rand(num_dof(sys),1)*
					0.01))
				sol_defl, = newton(F,J,ic_main_p,bias,
					optnewton_defl,deflation_op)
				vr_branch = vr_deflation .+ vr_main[i]
				bias_branch = hcat(zeros(length(vr_branch)),
					vr_branch)
				ic_branch = reshape(sol_defl,sys)
				j_defl, = non_equilib(d,bias_branch,ic_branch,
					int_contacts=2)
				plot(vr_branch,j_defl)
				ylim(0,1e-9)
				push!(deflation_op,sol_defl)
			catch e
				println("No branch found")
				break
			end
		end

	end

	return (nothing)

end

end

