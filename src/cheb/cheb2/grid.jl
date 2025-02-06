"""
    ChebyshevLobattoGrid{TF} <: AbstractGrid{TF}

The extrema of Chebyshev polynomials are called the Chebyshev points of the second kind, or Chebyshev extreme points, or Chebyshev–Lobatto points.
"""
struct ChebyshevLobattoGrid{TF<:AbstractFloat} <: AbstractGrid{TF}
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

export ChebyshevLobattoGrid
