abstract type AbstractChebyshevNode end

struct ChebyshevFirstKindNode <: AbstractChebyshevNode end
struct ChebyshevSecondKindNode <: AbstractChebyshevNode end

"""
    ChebyshevGrid{TF,TNode}(n, lower_bound, upper_bound, kind=2) where TF<:AbstractFloat

Build a Chebyshev grid of size `n` in the interval `[lower_bound, upper_bound]`. The grid can be of the first or second kind.

# Arguments
- `n::Integer`: Number of grid points
- `lower_bound::TF`: Lower bound of the interval
- `upper_bound::TF`: Upper bound of the interval
- `kind::Integer`: Kind of Chebyshev grid (1 or 2)
"""
struct ChebyshevGrid{TF<:AbstractFloat,TNode<:AbstractChebyshevNode} <: AbstractGrid{TF}
    lower_bound::TF  # Lower bound of the interval
    upper_bound::TF  # Upper bound of the interval
    data::Vector{TF} # Grid points
    node::TNode      # Node type

    function ChebyshevGrid(
        n::Integer, lower_bound::TF, upper_bound::TF; kind::Integer=2
    ) where {TF<:AbstractFloat}
        @argcheck n >= 0 "n must be nonnegative"
        @argcheck upper_bound > lower_bound "upper_bound must be greater than lower_bound"
        @argcheck kind == 1 || kind == 2 "kind must be either 1 or 2"

        if kind == 1
            TNode = ChebyshevFirstKindNode
        elseif kind == 2
            TNode = ChebyshevSecondKindNode
        end

        node = TNode()
        data = _cheb_points(node, TF, n, lower_bound, upper_bound)
        return new{TF,TNode}(lower_bound, upper_bound, data, node)
    end
end

Base.length(grid::ChebyshevGrid) = length(grid.data)
Base.size(grid::ChebyshevGrid) = size(grid.data)

Base.@propagate_inbounds function Base.getindex(grid::ChebyshevGrid, i)
    return grid.data[i]
end

Base.@propagate_inbounds function Base.iterate(grid::ChebyshevGrid, state...)
    return iterate(grid.data, state...)
end

Base.keys(grid::ChebyshevGrid) = keys(grid.data)
Base.firstindex(grid::ChebyshevGrid) = firstindex(grid.data)
Base.lastindex(grid::ChebyshevGrid) = lastindex(grid.data)
Base.eltype(::Type{ChebyshevGrid{TF,TNode}}) where {TF,TNode} = TF
Base.IndexStyle(::Type{ChebyshevGrid}) = IndexLinear()

export ChebyshevGrid
