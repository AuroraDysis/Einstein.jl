module FDMSuite

using ..Utils

using ArgCheck: @argcheck

using FillArrays
using FastBroadcast
using LinearAlgebra
using BandedMatrices

include("utils.jl")
include("grid.jl")
include("weights.jl")
include("stencil.jl")
include("diss.jl")
include("fdmop.jl")
include("diffmat.jl")
include("dissmat.jl")
include("integrate.jl")
# include("interp.jl")

end
