abstract type AbstractChebyshevNode end

"""
    ChebyshevFirstKindNode <: AbstractChebyshevNode

The Chebyshev nodes of the first kind, also called the Chebyshev zeros, are the zeros of the Chebyshev polynomials of the first kind.
"""
struct ChebyshevFirstKindNode <: AbstractChebyshevNode end

"""
    ChebyshevSecondKindNode <: AbstractChebyshevNode

The Chebyshev nodes of the second kind are also referred to as Chebyshev–Lobatto points or Chebyshev extreme points.
"""
struct ChebyshevSecondKindNode <: AbstractChebyshevNode end

struct ChebyshevGrid{TF<:AbstractFloat,TNode<:AbstractChebyshevNode} <: AbstractGrid{TF}
    lower_bound::TF  # Lower bound of the interval
    upper_bound::TF  # Upper bound of the interval
    data::Vector{TF} # Grid points
    type::TNode      # Type of Chebyshev nodes
end

"""
    cheb_grid([node=ChebyshevSecondKindNode()], n, lower_bound, upper_bound)

Build a Chebyshev grid of size `n` in the interval `[lower_bound, upper_bound]`.

# Arguments
- `TNode`: Type of Chebyshev nodes, either `ChebyshevFirstKindNode` or `ChebyshevSecondKindNode`
- `n::Integer`: Number of grid points
- `lower_bound::TF`: Lower bound of the interval
- `upper_bound::TF`: Upper bound of the interval
"""
function cheb_grid(
    type::ChebyshevFirstKindNode, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    @argcheck n >= 0 "n must be nonnegative"
    @argcheck upper_bound > lower_bound "upper_bound must be greater than lower_bound"
    data = cheb1_pts(TF, n, lower_bound, upper_bound)
    return ChebyshevGrid{TF,ChebyshevFirstKindNode}(lower_bound, upper_bound, data, type)
end

function cheb_grid(
    type::ChebyshevSecondKindNode, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    @argcheck n >= 0 "n must be nonnegative"
    @argcheck upper_bound > lower_bound "upper_bound must be greater than lower_bound"
    data = cheb2_pts(TF, n, lower_bound, upper_bound)
    return ChebyshevGrid{TF,ChebyshevSecondKindNode}(lower_bound, upper_bound, data, type)
end

function cheb_grid(n::Integer, lower_bound::TF, upper_bound::TF) where {TF<:AbstractFloat}
    return cheb_grid(ChebyshevSecondKindNode(), n, lower_bound, upper_bound)
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

export ChebyshevGrid, cheb_grid, ChebyshevFirstKindNode, ChebyshevSecondKindNode
