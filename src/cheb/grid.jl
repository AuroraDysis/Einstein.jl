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

@enumx ChebyshevNode FirstKind SecondKind

struct ChebyshevGrid{TF<:AbstractFloat,TNode<:AbstractChebyshevNode} <: AbstractGrid{TF}
    lower_bound::TF  # Lower bound of the interval
    upper_bound::TF  # Upper bound of the interval
    data::Vector{TF} # Grid points
end

"""
    cheb_grid(n::Integer, lower_bound::TF, upper_bound::TF; node::ChebyshevNode.T=ChebyshevNode.SecondKind)

Build a Chebyshev grid of size `n` in the interval `[lower_bound, upper_bound]`.

# Arguments
- `n::Integer`: Number of grid points
- `lower_bound::TF`: Lower bound of the interval
- `upper_bound::TF`: Upper bound of the interval
- `node::ChebyshevNode.T`: Specifies the Chebyshev node type (FirstKind or SecondKind)
"""
function cheb_grid(
    n::Integer,
    lower_bound::TF,
    upper_bound::TF;
    node::ChebyshevNode.T=ChebyshevNode.SecondKind,
) where {TF<:AbstractFloat}
    @argcheck n >= 0 "n must be nonnegative"
    @argcheck upper_bound > lower_bound "upper_bound must be greater than lower_bound"

    if node == ChebyshevNode.FirstKind
        TNode = ChebyshevFirstKindNode
        data = cheb1_pts(TF, n, lower_bound, upper_bound)
    elseif node == ChebyshevNode.SecondKind
        TNode = ChebyshevSecondKindNode
        data = cheb2_pts(TF, n, lower_bound, upper_bound)
    else
        throw(ArgumentError("Invalid Chebyshev node type: $node"))
    end

    return ChebyshevGrid{TF,TNode}(lower_bound, upper_bound, data)
end

# Overloads for ChebyshevGrid
Base.length(grid::ChebyshevGrid) = length(grid.data)
Base.size(grid::ChebyshevGrid) = size(grid.data)

Base.@propagate_inbounds function Base.getindex(grid::ChebyshevGrid, i)
    return grid.data[i]
end

# Additional Base overload: iteration
function Base.iterate(grid::ChebyshevGrid, state=(eachindex(grid.data),))
    return iterate(grid.data, state)
end

Base.keys(grid::ChebyshevGrid) = keys(grid.data)
Base.firstindex(grid::ChebyshevGrid) = firstindex(grid.data)
Base.lastindex(grid::ChebyshevGrid) = lastindex(grid.data)
Base.eltype(::Type{ChebyshevGrid{TF}}) where {TF<:AbstractFloat} = TF
Base.IndexStyle(::Type{ChebyshevGrid}) = IndexLinear()

export ChebyshevNode, ChebyshevGrid, cheb_grid
