"""
    ChebyshevInterpolation(grid::ChebyshevGrid{TF,Basis}) where {TF<:AbstractFloat,Basis<:AbstractBasisFunction}
    (itp::ChebyshevInterpolation{TF})(values::AbstractVector{TFC}, x::TF) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}

Construct a barycentric interpolation with precomputed weights for a Chebyshev grid.
"""
struct ChebyshevInterpolation{TF<:AbstractFloat,Basis<:AbstractBasisFunction}
    grid::ChebyshevGrid{TF,Basis}
    weights::Vector{TF}

    function ChebyshevInterpolation(
        grid::ChebyshevGrid{TF,Basis}
    ) where {TF<:AbstractFloat,Basis<:AbstractBasisFunction}
        weights = _cheb_barycentric_weights(grid.basis, TF, length(grid))
        return new{TF,Basis}(grid, weights)
    end
end

function (itp::ChebyshevInterpolation{TF})(
    values::AbstractVector{TFC}, x::TF
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    (; grid, weights) = itp
    @argcheck grid.lower_bound <= x <= grid.upper_bound "x is out of range"

    points = grid.data
    return barycentric_interpolate(x, points, values, weights)
end

export ChebyshevInterpolation
