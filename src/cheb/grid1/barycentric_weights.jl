"""
    chebtech1_barycentric_weights([TF=Float64], n::Integer) where {TF<:AbstractFloat}

Compute the barycentric weights for Chebyshev points of the 1st kind.

# Arguments
- `TF`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# References

- [berrut2004barycentric](@citet*)
- [chebfun/@chebtech1/barywts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech1/barywts.m)

See also: [`BarycentricInterpolation`](@ref), [`chebtech1_points`](@ref)
"""
function chebtech1_barycentric_weights(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    if n == 0
        return TF[]
    elseif n == 1
        return TF[one(TF)]
    end

    half = one(TF) / 2
    pi_over_n = convert(TF, π) / n
    v = Vector{TF}(undef, n)
    @inbounds for j in 0:(n - 1)
        θ = (n - j - half) * pi_over_n
        v[j + 1] = sin(θ)
    end

    # The following flipping trick forces symmetry. Also due to the nature of 
    # the sine function, those computed with a big argument are replaced by ones
    # with a small argument, improving the relative accuracy.
    half_n = floor(typeof(n), n / 2)
    # Copy values from end to beginning for symmetry
    @inbounds for i in 1:half_n
        v[i] = v[n - i + 1]
    end

    # Flip signs for odd indices
    @inbounds for i in (n - 1):-2:1
        v[i] = -v[i]
    end

    return v
end

function chebtech1_barycentric_weights(n::Integer)
    return chebtech1_barycentric_weights(Float64, n)
end
