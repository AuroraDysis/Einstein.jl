@doc raw"""
    barycentric_weights(x::AbstractVector{TF}) where {TF<:AbstractFloat}
    barycentric_weights(x::AbstractRange{TF}) where {TF<:AbstractFloat}
    barycentric_weights(::Type{TF}, n::Integer) where {TF<:AbstractFloat}

Computes barycentric weights for the nodes `x`.
Alternatively, for $n + 1$ equidistant nodes, the weights can be computed by passing `n`.

# References
- [Berrut2004](@citet*)
"""
function barycentric_weights(x::AbstractVector{TF}) where {TF<:AbstractFloat}
    return TF[inv(prod(x[i] - x[j] for j in eachindex(x) if j != i)) for i in eachindex(x)]
end

function barycentric_weights(x::AbstractRange{TF}) where {TF<:AbstractFloat}
    n = length(x) - 1
    return barycentric_weights(TF, n)
end

function barycentric_weights(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    return TF[(-1)^j * binomial(n, j) for j in 0:n]
end

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
struct BarycentricInterpolation{TF<:AbstractFloat}
    points::AbstractVector{TF}   # Grid points
    weights::AbstractVector{TF}  # Barycentric weights

    function BarycentricInterpolation(points::AbstractVector{TF}) where {TF<:AbstractFloat}
        weights = barycentric_weights(points)
        return new{TF}(points, weights)
    end
end

function (itp::BarycentricInterpolation{TF})(
    values::AbstractVector{TFC}, x::TF
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    (; points, weights) = itp
    @argcheck points[begin] <= x <= points[end] "x is out of range"

    return barycentric_interpolate(x, points, values, weights)
end

export BarycentricInterpolation, barycentric_weights, barycentric_interpolate
