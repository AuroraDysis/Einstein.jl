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

"""
    ChebyshevGrid{TF <: AbstractFloat} <: AbstractGrid{TF}

The zeros of Chebyshev polynomials are called Chebyshev points of the first kind, Chebyshev nodes, or, more formally, Chebyshev–Gauss points.
"""
struct ChebyshevGrid{TF<:AbstractFloat,TNode<:AbstractChebyshevNode} <: AbstractGrid{TF}
    x_min::TF # Lower bound of the interval
    x_max::TF # Upper bound of the interval
    n::Integer # Number of grid points
    x::Vector{TF} # Grid points
end

function cheb_grid(
    ::Type{TNode}, x_min::TF, x_max::TF, n::Integer
) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    @argcheck n >= 0 "n must be nonnegative"
    @argcheck x_max > x_min "x_max must be greater than x_min"
    return ChebyshevGrid{TF,TNode}(; x_min=x_min, x_max=x_max, n=n)
end

Base.length(grid::ChebyshevGrid) = grid.n
Base.size(grid::ChebyshevGrid) = (grid.n,)
Base.@propagate_inbounds Base.getindex(grid::ChebyshevGrid, i::Integer) = grid.x[i]

export ChebyshevGrid
