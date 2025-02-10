"""
    ChebyshevNode

FirstKind=1, The Chebyshev nodes of the first kind, also called the Chebyshev zeros, are the zeros of the Chebyshev polynomials of the first kind.
SecondKind=2, The Chebyshev nodes of the second kind are also referred to as Chebyshev–Lobatto points or Chebyshev extreme points.
"""
@enumx ChebyshevNode FirstKind = 1 SecondKind = 2

struct ChebyshevGrid{TF<:AbstractFloat} <: AbstractGrid{TF}
    lower_bound::TF  # Lower bound of the interval
    upper_bound::TF  # Upper bound of the interval
    data::Vector{TF} # Grid points
    type::ChebyshevNode.T
end

"""
    cheb_grid(n, lower_bound, upper_bound, kind=ChebyshevNode.SecondKind)

Build a Chebyshev grid of size `n` in the interval `[lower_bound, upper_bound]`.

# Arguments
- `n::Integer`: Number of grid points
- `lower_bound::TF`: Lower bound of the interval
- `upper_bound::TF`: Upper bound of the interval
"""
function cheb_grid(
    n::Integer,
    lower_bound::TF,
    upper_bound::TF,
    kind::ChebyshevNode.T=ChebyshevNode.SecondKind,
) where {TF<:AbstractFloat}
    @argcheck n >= 0 "n must be nonnegative"
    @argcheck upper_bound > lower_bound "upper_bound must be greater than lower_bound"
    if kind == ChebyshevNode.FirstKind
        data = cheb1_pts(TF, n, lower_bound, upper_bound)
    elseif kind == ChebyshevNode.SecondKind
        data = cheb2_pts(TF, n, lower_bound, upper_bound)
    else
        throw(
            ArgumentError(
                "kind must be either ChebyshevNode.FirstKind or ChebyshevNode.SecondKind"
            ),
        )
    end
    return ChebyshevGrid{TF}(lower_bound, upper_bound, data, kind)
end

# Overloads for ChebyshevGrid
Base.length(grid::ChebyshevGrid) = length(grid.data)
Base.size(grid::ChebyshevGrid) = size(grid.data)

Base.@propagate_inbounds function Base.getindex(grid::ChebyshevGrid, i)
    return grid.data[i]
end

# Additional Base overload: iteration
Base.@propagate_inbounds function Base.iterate(grid::ChebyshevGrid, state...)
    return iterate(grid.data, state...)
end

Base.keys(grid::ChebyshevGrid) = keys(grid.data)
Base.firstindex(grid::ChebyshevGrid) = firstindex(grid.data)
Base.lastindex(grid::ChebyshevGrid) = lastindex(grid.data)
Base.eltype(::Type{ChebyshevGrid{TF}}) where {TF<:AbstractFloat} = TF
Base.IndexStyle(::Type{ChebyshevGrid}) = IndexLinear()

export ChebyshevGrid, cheb_grid, ChebyshevNode
