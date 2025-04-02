"""
    cheb_lobatto_barycentric_weights([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Compute the barycentric weights for Chebyshev points of the 2nd kind.

# Arguments
- `TF`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# References

- [chebfun/@chebtech2/barywts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/barywts.m)

See also: [`BarycentricInterpolation`](@ref)
"""
function cheb_lobatto_barycentric_weights(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    if n == 0
        return TF[]
    elseif n == 1
        return ones(TF, 1)
    end

    weights = ones(TF, n)

    @inbounds begin
        half = one(TF) / 2
        weights[(end - 1):-2:1] .= -1
        weights[1] *= half
        weights[end] = half
    end

    return weights
end

function cheb_lobatto_barycentric_weights(n::Integer)
    return cheb_lobatto_barycentric_weights(Float64, n)
end

export cheb_lobatto_barycentric_weights
