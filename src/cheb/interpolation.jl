"""
    ChebyshevInterpolation(grid::ChebyshevGrid{TF,TNode}) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    (itp::ChebyshevInterpolation{TF})(values::AbstractVector{TFC}, x::TF) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}

Construct a barycentric interpolation with precomputed weights for a Chebyshev grid.
"""
struct ChebyshevInterpolation{TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    grid::ChebyshevGrid{TF,TNode}
    weights::Vector{TF}

    function ChebyshevInterpolation(
        grid::ChebyshevGrid{TF,TNode}
    ) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
        weights = _cheb_barycentric_weights(grid.node, TF, length(grid))
        return new{TF,TNode}(grid, weights)
    end
end

function (itp::ChebyshevInterpolation{TF})(
    values::AbstractVector{TFC}, x::TF
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    (; grid, weights) = itp
    points = grid.data
    return barycentric_interpolate(x, points, values, weights)
end

export ChebyshevInterpolation
