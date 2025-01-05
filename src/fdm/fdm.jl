module FDMSuite

using ..Utils

using ArgCheck: @argcheck

using FillArrays
using LinearAlgebra
using BandedMatrices

include("utils.jl")
include("grid.jl")
include("weights.jl")
include("stencil.jl")
include("diss.jl")
include("fdmop.jl")
include("diffmat.jl")

end
