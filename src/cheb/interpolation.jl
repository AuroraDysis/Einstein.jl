"""
    ChebyshevInterpolation(grid::ChebyshevGrid{TF,Basis}) where {TF<:AbstractFloat,Basis<:AbstractBasisFunction}
    (itp::ChebyshevInterpolation{TF})(values::AbstractVector{TFC}, x::TF) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}

Construct a barycentric interpolation with precomputed weights for a Chebyshev grid.
"""
struct ChebyshevInterpolation{TF<:AbstractFloat}
    points::Vector{TF}
    weights::Vector{TF}

    function ChebyshevInterpolation(
        points::AbstractVector{TF}, weights::AbstractVector{TFC}
    ) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
        return new{TF}(points, weights)
    end
end

function (itp::ChebyshevInterpolation{TF})(
    values::AbstractVector{TFC}, x::TF
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    (; points, weights) = itp

    @boundscheck begin
        @argcheck minimum(points) <= x <= maximum(points) "x is out of range"
    end

    return barycentric_interpolate(x, points, values, weights)
end

export ChebyshevInterpolation
