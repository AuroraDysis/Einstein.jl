"""
    barycentric_weights([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Compute the barycentric weights for Chebyshev points of the 2nd kind.

# Arguments
- `TF`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# References

- [chebfun/@chebtech2/barywts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/barywts.m)

See also: [`BarycentricInterpolation`](@ref), [`points`](@ref)
"""
function barycentric_weights(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    if n == 0
        return TF[]
    elseif n == 1
        return TF[one(TF)]
    end

    w = ones(TF, n)

    @inbounds begin
        half = one(TF) / 2
        w[(end - 1):-2:1] .= -1
        w[1] *= half
        w[end] = half
    end

    return w
end

function barycentric_weights(n::Integer)
    return barycentric_weights(Float64, n)
end

export barycentric_weights
