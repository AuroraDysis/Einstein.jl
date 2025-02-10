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

module ChebyshevNode
using ..ChebSuite: ChebyshevFirstKindNode, ChebyshevSecondKindNode

const FirstKind = ChebyshevFirstKindNode
const SecondKind = ChebyshevSecondKindNode
end

struct ChebyshevGrid{TF<:AbstractFloat,TNode<:AbstractChebyshevNode} <: AbstractGrid{TF}
    lower_bound::TF  # Lower bound of the interval
    upper_bound::TF  # Upper bound of the interval
    data::Vector{TF} # Grid points
end

"""
    cheb_grid([TNode=ChebyshevNode.SecondKind], n, lower_bound, upper_bound)

Build a Chebyshev grid of size `n` in the interval `[lower_bound, upper_bound]`.

# Arguments
- `TNode`: Type of Chebyshev nodes, either `ChebyshevNode.FirstKind` or `ChebyshevNode.SecondKind`
- `n::Integer`: Number of grid points
- `lower_bound::TF`: Lower bound of the interval
- `upper_bound::TF`: Upper bound of the interval
"""
function cheb_grid(
    ::Type{TNode}, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    @argcheck n >= 0 "n must be nonnegative"
    @argcheck upper_bound > lower_bound "upper_bound must be greater than lower_bound"

    if TNode == ChebyshevNode.FirstKind
        data = cheb1_pts(TF, n, lower_bound, upper_bound)
    elseif TNode == ChebyshevNode.SecondKind
        data = cheb2_pts(TF, n, lower_bound, upper_bound)
    else
        throw(ArgumentError("Invalid Chebyshev node type: $TNode"))
    end

    return ChebyshevGrid{TF,TNode}(lower_bound, upper_bound, data)
end

function cheb_grid(n::Integer, lower_bound::TF, upper_bound::TF) where {TF<:AbstractFloat}
    return cheb_grid(ChebyshevNode.SecondKind, n, lower_bound, upper_bound)
end

# Overloads for ChebyshevGrid
Base.length(grid::ChebyshevGrid) = length(grid.data)
Base.size(grid::ChebyshevGrid) = size(grid.data)

Base.@propagate_inbounds function Base.getindex(grid::ChebyshevGrid, i)
    return grid.data[i]
end

# Additional Base overload: iteration
function Base.iterate(grid::ChebyshevGrid, i=firstindex(grid.data))
    return iterate(grid.data, i)
end

Base.keys(grid::ChebyshevGrid) = keys(grid.data)
Base.firstindex(grid::ChebyshevGrid) = firstindex(grid.data)
Base.lastindex(grid::ChebyshevGrid) = lastindex(grid.data)
Base.eltype(::Type{ChebyshevGrid{TF}}) where {TF<:AbstractFloat} = TF
Base.IndexStyle(::Type{ChebyshevGrid}) = IndexLinear()

export ChebyshevNode, ChebyshevGrid, cheb_grid
