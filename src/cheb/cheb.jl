module ChebSuite

import AbstractFFTs: Plan
using ..Utils
using FFTW
using FastTransforms
using ArgCheck: @argcheck
using FastBroadcast: @..
using LinearAlgebra
using FillArrays
using SparseArrays
using BandedMatrices
using ToeplitzMatrices: Toeplitz, Hankel

"""
    ChebyshevGaussGrid{TF} <: AbstractGrid{TF}

The zeros of Chebyshev polynomials are called Chebyshev points of the first kind, Chebyshev nodes, or, more formally, Chebyshev–Gauss points.
"""
struct ChebyshevGaussGrid{TF} <: AbstractGrid{TF} where {TF<:AbstractFloat}
    x_min::TF # Lower bound of the interval
    x_max::TF # Upper bound of the interval
    n::Integer # Number of grid points
    grid::Vector{TF} # Grid points

    function ChebyshevGaussGrid(x_min::TF, x_max::TF, n::Integer) where {TF<:AbstractFloat}
        @argcheck n >= 0 "n must be nonnegative"
        @argcheck x_max > x_min "x_max must be greater than x_min"

        grid = cheb1_pts(TF, n, x_min, x_max)
        return new{TF}(x_min, x_max, n, grid)
    end
end

Base.length(grid::ChebyshevGaussGrid) = grid.n
Base.size(grid::ChebyshevGaussGrid) = (grid.n,)
@propagate_inbounds Base.getindex(grid::ChebyshevGaussGrid, i::Integer) = grid.grid[i]

"""
    ChebyshevLobattoGrid{TF} <: AbstractGrid{TF}

The extrema of Chebyshev polynomials are called the Chebyshev points of the second kind, or Chebyshev extreme points, or Chebyshev–Lobatto points.
"""
struct ChebyshevLobattoGrid{TF} <: AbstractGrid{TF} where {TF<:AbstractFloat}
    x_min::TF # Lower bound of the interval
    x_max::TF # Upper bound of the interval
    n::Integer # Number of grid points
    grid::Vector{TF} # Grid points

    function ChebyshevLobattoGrid(
        x_min::TF, x_max::TF, n::Integer
    ) where {TF<:AbstractFloat}
        @argcheck n >= 0 "n must be nonnegative"
        @argcheck x_max > x_min "x_max must be greater than x_min"

        grid = cheb2_pts(TF, n, x_min, x_max)
        return new{TF}(x_min, x_max, n, grid)
    end
end

Base.length(grid::ChebyshevLobattoGrid) = grid.n
Base.size(grid::ChebyshevLobattoGrid) = (grid.n,)
@propagate_inbounds Base.getindex(grid::ChebyshevLobattoGrid, i::Integer) = grid.grid[i]

export ChebyshevGaussGrid, ChebyshevLobattoGrid

include("utils/bary.jl")
include("utils/clenshaw.jl")
include("utils/feval.jl")
include("utils/fft.jl")
include("utils/diff.jl")
include("utils/cumsum.jl")
include("utils/bary_diffmat.jl")
include("utils/cumsummat.jl")
include("utils/dissipation.jl")

include("cheb1/angles.jl")
include("cheb1/pts.jl")
include("cheb1/coeffs2vals.jl")
include("cheb1/vals2coeffs.jl")
include("cheb1/barywts.jl")
include("cheb1/interp.jl")
include("cheb1/quad.jl")
include("cheb1/diffmat.jl")
include("cheb1/asmat.jl")
include("cheb1/cumsummat.jl")

include("cheb2/angles.jl")
include("cheb2/pts.jl")
include("cheb2/coeffs2vals.jl")
include("cheb2/vals2coeffs.jl")
include("cheb2/barywts.jl")
include("cheb2/interp.jl")
include("cheb2/quad.jl")
include("cheb2/diffmat.jl")
include("cheb2/asmat.jl")
include("cheb2/cumsummat.jl")

# Rectangular spectral collocation
include("rect/rectint.jl")
include("rect/rectdiff.jl")

# Ultraspherical spectral method
include("ultra/diffmat.jl")
include("ultra/convertmat.jl")
include("ultra/sphankel.jl")
include("ultra/multmat.jl")

end
