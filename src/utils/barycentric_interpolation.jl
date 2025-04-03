@doc raw"""
    barycentric_interpolate(x::TF, points::AbstractVector{TF}, values::AbstractVector{TFC}, weights::Vector{TF}) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}

Evaluate the value of the barycentric interpolation formula with nodes `points`, function values `values`
and barycentric weights `weights` at point `x`.

# References
- [chebfun/@chebtech2/bary.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/bary.m)
- [Berrut2004](@citet*)
"""
function barycentric_interpolate(
    x::TF, points::AbstractVector{TF}, values::AbstractVector{TFC}, weights::Vector{TF}
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    @boundscheck begin
        @argcheck length(points) > 1 "points must have at least two elements"
        @argcheck length(points) == length(values) "points and values must have the same length"
        @argcheck length(points) == length(weights) "points and weights must have the same length"
    end

    p = zero(TFC)
    q = zero(TFC)

    @inbounds for i in eachindex(points)
        Δx = x - points[i]

        if iszero(Δx)
            return values[i]
        end

        wi = weights[i] / Δx
        p += wi * values[i]
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
    (itp::BarycentricInterpolation{TF})(values::AbstractVector{TFC}, x::TF) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}

Evaluate the interpolant at point `x` for function values.

# Reference
- [chebfun/@chebtech2/bary.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/bary.m)
- [Berrut2004](@citet*)
"""
struct BarycentricInterpolation{
    TF<:AbstractFloat,VT1<:AbstractVector{TF},VT2<:AbstractVector{TF}
}
    points::VT1   # Grid points
    weights::VT2  # Barycentric weights

    function BarycentricInterpolation(
        points::VT1, weights::VT2
    ) where {TF<:AbstractFloat,VT1<:AbstractVector{TF},VT2<:AbstractVector{TF}}
        return new{TF,VT1,VT2}(points, weights)
    end
end

function (itp::BarycentricInterpolation{TF,VT1,VT2})(
    values::AbstractVector{TFC}, x::TF
) where {
    TF<:AbstractFloat,
    TFC<:Union{TF,Complex{TF}},
    VT1<:AbstractVector{TF},
    VT2<:AbstractVector{TF},
}
    (; points, weights) = itp

    @boundscheck begin
        @argcheck first(points) <= x <= last(points) "x is out of range"
    end

    return barycentric_interpolate(x, points, values, weights)
end

export BarycentricInterpolation, barycentric_interpolate
