function barycentric_kernal(
    x::TF, points::AbstractVector{TF}, ys::AbstractVector{TR}, weights::Vector{TF}
) where {TF<:AbstractFloat,TR<:Union{TF,Complex{TF}}}
    p = zero(TR)
    q = zero(TR)

    @inbounds for i in eachindex(points)
        Δx = x - points[i]

        if iszero(Δx)
            return ys[i]
        end

        wi = weights[i] / Δx
        p += wi * ys[i]
        q += wi
    end

    return p / q
end

"""
    BarycentricInterpolation{TF<:AbstractFloat}(points::AbstractVector{TF}, weights::AbstractVector{TF})

A structure representing barycentric interpolation with precomputed weights.

# Fields
- `points::AbstractVector{TF}`: Vector of interpolation points (typically Chebyshev points)
- `weights::AbstractVector{TF}`: Vector of barycentric weights

# Methods
    (itp::BarycentricInterpolation{TF})(values::AbstractVector{TR}, x::TF) where {TF<:AbstractFloat,TR<:Union{TF,Complex{TF}}}

Evaluate the interpolant at point `x` for function values.

# Reference
- [chebfun/@chebtech2/bary.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/bary.m)
"""
struct BarycentricInterpolation{TF<:AbstractFloat}
    points::AbstractVector{TF}   # Grid points
    weights::AbstractVector{TF}        # Barycentric weights
end

function (itp::BarycentricInterpolation{TF})(
    values::AbstractVector{TR}, x::TF
) where {TF<:AbstractFloat,TR<:Union{TF,Complex{TF}}}
    (; points, weights) = itp

    return barycentric_kernal(x, points, values, weights)
end

export BarycentricInterpolation
