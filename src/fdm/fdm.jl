module FDMSuite

using ..Utils

using ArgCheck: @argcheck

using FillArrays
using FastBroadcast
using LinearAlgebra
using BandedMatrices

"""
    UniformGrid{TF} <: AbstractGrid{TF}

Uniform grid with constant spacing.
"""
struct UniformGrid{TF} <: AbstractGrid{TF} where {TF<:AbstractFloat}
    x_min::TF # lower bound of the grid
    x_max::TF # upper bound of the grid
    Δx::TF # Grid spacing
    n::Integer # Number of grid points
    grid::StepRangeLen{TF}

    function UniformGrid(x_min::TF, x_max::TF, n::Integer) where {TF<:AbstractFloat}
        @argcheck n >= 0 "n must be nonnegative"
        @argcheck x_max > x_min "x_max must be greater than x_min"

        grid = range(x_min; stop=x_max, length=n)
        dx = step(grid)
        n = length(grid)

        return new{TF}(x_min, x_max, dx, n, grid)
    end
end

export UniformGrid

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
