"""
    ChebyshevInterpolation(grid::ChebyshevGrid{TF,TNode}) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}

Construct a barycentric interpolation with precomputed weights for a Chebyshev grid.
"""
struct ChebyshevInterpolation{TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    grid::ChebyshevGrid{TF,TNode}
    bary_itp::BarycentricInterpolation{TF}

    function ChebyshevInterpolation(
        grid::ChebyshevGrid{TF,TNode}
    ) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
        weights = _cheb_barycentric_weights(grid.node, TF, length(grid))
        return new{TF,TNode}(grid, BarycentricInterpolation(grid.data, weights))
    end
end

function (itp::ChebyshevInterpolation{TF})(
    f::AbstractVector{TFC}, x0::TF
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    return itp.bary_itp(f, x0)
end

export ChebyshevInterpolation
