@doc raw"""
    barycentric_weights(x::AbstractVector{TF}) where {TF<:AbstractFloat}
    barycentric_weights(x::AbstractRange{TF}) where {TF<:AbstractFloat}
    barycentric_weights(::Type{TF}, n::Integer) where {TF<:AbstractFloat}

Computes barycentric weights for the nodes `x`.
Alternatively, for $n + 1$ equidistant nodes, the weights can be computed by passing `n`.

The weights are scaled such that their infinity norm equals 1.

# References
- [Berrut2004](@citet*)
"""
function barycentric_weights(x::AbstractVector{TF}) where {TF<:AbstractFloat}
    @boundscheck begin
        @argcheck length(x) > 1 "x must have at least two elements"
    end

    weights = similar(x)

    x_min, x_max = extrema(x)
    if x_min ≈ x_max
        @warn "Input points are not distinct. Returning NaN weights."
        fill!(weights, NaN)
        return weights
    end

    # Capacity of the interval - helps prevent underflow/overflow
    c = 4 / (x_max - x_min)

    diff = similar(x)
    for i in eachindex(x)
        @.. diff = c * (x[i] - x)
        diff[i] = 1

        # Log-Sum-Exp trick to avoid underflow/overflow
        prod_sign = prod(sign.(diff))
        @.. diff = log(abs(diff))
        weights[i] = 1 / (prod_sign * exp(sum(diff)))
    end

    weights ./= norm(weights, Inf)

    return weights
end

function barycentric_weights(x::AbstractRange{TF}) where {TF<:AbstractFloat}
    n = length(x) - 1
    return barycentric_weights(TF, n)
end

function barycentric_weights(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    @boundscheck begin
        @argcheck n > 0 "n must be a positive integer"
    end

    norm_factor = binomial(n, n ÷ 2)
    weights = TF[(-1)^j * binomial(n, j)//norm_factor for j in 0:n]

    return weights
end

export barycentric_weights
