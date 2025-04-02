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
    @boundscheck begin
        @argcheck length(x) > 1 "x must have at least two elements"
    end

    return TF[inv(prod(x[i] - x[j] for j in eachindex(x) if j != i)) for i in eachindex(x)]
end

function barycentric_weights(x::AbstractRange{TF}) where {TF<:AbstractFloat}
    n = length(x) - 1
    return barycentric_weights(TF, n)
end

function barycentric_weights(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    @boundscheck begin
        @argcheck n > 0 "n must be a positive integer"
    end

    return TF[(-1)^j * binomial(n, j) for j in 0:n]
end

export barycentric_weights
