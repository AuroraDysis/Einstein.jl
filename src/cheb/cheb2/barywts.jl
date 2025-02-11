"""
    cheb2_barycentric_weights([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Compute the barycentric weights for Chebyshev points of the 2nd kind.

# Arguments
- `TF`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# References

- [chebfun/@chebtech2/barywts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/barywts.m)

See also: [`BarycentricInterpolation`](@ref), [`cheb2_points`](@ref)
"""
function cheb2_barycentric_weights(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
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

function cheb2_barycentric_weights(n::Integer)
    return cheb2_barycentric_weights(Float64, n)
end

function cheb_barycentric_weights(
    ::ChebyshevSecondKindNode, ::Type{TF}, n::Integer
) where {TF<:AbstractFloat}
    return cheb2_barycentric_weights(TF, n)
end

export cheb2_barycentric_weights
