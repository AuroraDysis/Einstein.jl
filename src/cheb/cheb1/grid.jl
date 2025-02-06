"""
    ChebyshevGaussGrid{TF <: AbstractFloat} <: AbstractGrid{TF}

The zeros of Chebyshev polynomials are called Chebyshev points of the first kind, Chebyshev nodes, or, more formally, Chebyshev–Gauss points.
"""
struct ChebyshevGaussGrid{TF<:AbstractFloat} <: AbstractGrid{TF}
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

export ChebyshevGaussGrid
