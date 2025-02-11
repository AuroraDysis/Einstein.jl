"""
    ChebyshevInterpolation(grid::ChebyshevGrid{TF}) where {TF<:AbstractFloat}

Construct a barycentric interpolation with precomputed weights for a Chebyshev grid.
"""
struct ChebyshevInterpolation{TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
    grid::ChebyshevGrid{TF,TNode}
    bary_itp::BarycentricInterpolation{TF}

    function ChebyshevInterpolation(
        grid::ChebyshevGrid{TF,TNode}
    ) where {TF<:AbstractFloat,TNode<:AbstractChebyshevNode}
        weights = cheb_barycentric_weights(grid.node, TF, length(grid))
        return new{TF,TNode}(grid, BarycentricInterpolation(grid.data, weights))
    end
end

function (itp::ChebyshevInterpolation{TF})(
    f::AbstractVector{TR}, x0::TF
) where {TF<:AbstractFloat,TR<:Union{TF,Complex{TF}}}
    return itp.bary_itp(f, x0)
end

export ChebyshevInterpolation
