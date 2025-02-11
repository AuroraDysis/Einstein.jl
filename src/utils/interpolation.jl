"""
    BarycentricInterpolation{TF<:AbstractFloat}(points::Vector{TF}, weights::Vector{TF})

A structure representing barycentric interpolation with precomputed weights.

# Fields
- `points::Vector{TF}`: Vector of interpolation points (typically Chebyshev points)
- `weights::Vector{TF}`: Vector of barycentric weights

# Methods
    (itp::BarycentricInterpolation{TF})(f::AbstractVector{TR}, x0::TF) where {TF<:AbstractFloat,TR<:Union{TF,Complex{TF}}}

Evaluate the interpolant at point `x0` for function values `f`.

# Reference
- [chebfun/@chebtech2/bary.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/bary.m)
"""
struct BarycentricInterpolation{TF<:AbstractFloat}
    points::Vector{TF} # Grid points
    weights::Vector{TF}  # Barycentric weights
end

function (itp::BarycentricInterpolation{TF})(
    f::AbstractVector{TR}, x0::TF
) where {TF<:AbstractFloat,TR<:Union{TF,Complex{TF}}}
    x = itp.points
    w = itp.weights

    p = zero(TR)
    q = zero(TR)

    @inbounds for i in eachindex(x)
        if x0 == x[i]
            return f[i]
        end

        wi = w[i] / (x0 - x[i])
        p += wi * f[i]
        q += wi
    end

    return p / q
end

export BarycentricInterpolation
