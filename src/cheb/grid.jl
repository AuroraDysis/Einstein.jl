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
    x_min::TF # Lower bound of the interval
    x_max::TF # Upper bound of the interval
    n::Integer # Number of grid points
    x::Vector{TF} # Grid points
end

function cheb_grid(
    x_min::TF, x_max::TF, n::Integer; node::ChebyshevNode.T=ChebyshevNode.SecondKind
) where {TF<:AbstractFloat}
    @argcheck n >= 0 "n must be nonnegative"
    @argcheck x_max > x_min "x_max must be greater than x_min"
    if node == ChebyshevNode.FirstKind
        TNode = ChebyshevFirstKindNode
        x = cheb1_pts(TF, n, x_min, x_max)
    elseif node == ChebyshevNode.SecondKind
        TNode = ChebyshevSecondKindNode
        x = cheb2_pts(TF, n, x_min, x_max)
    else
        throw(ArgumentError("Invalid Chebyshev node type: $node"))
    end

    return ChebyshevGrid{TF,TNode}(x_min, x_max, n, x)
end

Base.length(grid::ChebyshevGrid) = grid.n
Base.size(grid::ChebyshevGrid) = (grid.n,)
Base.@propagate_inbounds Base.getindex(grid::ChebyshevGrid, i::Integer) = grid.x[i]

export ChebyshevNode, ChebyshevGrid, cheb_grid
