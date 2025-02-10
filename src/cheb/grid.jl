abstract type AbstractChebyshevNode end

"""
    ChebyshevNodeFirstKind <: AbstractChebyshevNode

The Chebyshev nodes of the first kind, also called the Chebyshev zeros, are the zeros of the Chebyshev polynomials of the first kind.
"""
struct ChebyshevFirstKindNode <: AbstractChebyshevNode end

"""
    ChebyshevNodeSecondKind <: AbstractChebyshevNode

The Chebyshev nodes of the second kind are also referred to as Chebyshev-Lobatto points or Chebyshev extreme points.
"""
struct ChebyshevSecondKindNode <: AbstractChebyshevNode end

@enumx ChebyshevNode FirstKind SecondKind

struct ChebyshevGrid{TF<:AbstractFloat,TNode<:AbstractChebyshevNode} <: AbstractGrid{TF}
    lower_bound::TF # Lower bound of the interval
    upper_bound::TF # Upper bound of the interval
    data::Vector{TF} # Grid points
end

function cheb_grid(
    lower_bound::TF,
    upper_bound::TF,
    n::Integer;
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

Base.length(grid::ChebyshevGrid) = length(grid.data)
Base.size(grid::ChebyshevGrid) = size(grid.data)
Base.@propagate_inbounds Base.getindex(grid::ChebyshevGrid, i) = grid.data[i]

export ChebyshevNode, ChebyshevGrid, cheb_grid
