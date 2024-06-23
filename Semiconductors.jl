# Semiconductors module (skeleton)
# Sam Chinnery
# 2022-03-15

module Semiconductors

using Revise
using VoronoiFVM
using LinearAlgebra
using SparseArrays
using Symbolics
using PyPlot
using Printf

using GridVisualize
using ExtendableGrids
using ExtendableSparse
using SimplexGridFactory
using Triangulate

using ForwardDiff
using ReverseDiff
using Zygote

include("device_solvers_v13.jl")
export equilib
export non_equilib_sys
export non_equilib
export do_integrate

include("device_objects_v2.jl")

include("geometries_v10.jl")
export diode1d
export diode2d
export thy1d
export npn1
export npn2
export mos1

#include("iv_continuation_v5.jl")
#export IVState
#export IVCurve
#export initialize!
#export continuation!
#export current

include("vfvm_bk_interface_v1.jl")
export bk_residual
export bk_jacobian

include("current_bc_v2.jl")
export NewtonState
export NewtonParams
export current
export newton_current!
export verbose_noverbose
export verbose_print_iters

include("vfvm_flux_v4.jl")
export charge_neutrality_ic!
export newton!
export flux_integrate!

# module Semiconductors
end

