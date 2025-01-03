"""
    cheb2_barywts([T=Float64], n::Integer) where {T<:AbstractFloat}

Compute the barycentric weights for Chebyshev points of the 2nd kind.

# Arguments
- `T`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# References

- [chebfun/@chebtech2/barywts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/barywts.m)

See also: [`bary`](@ref), [`cheb2_pts`](@ref)
"""
function cheb2_barywts(::Type{T}, n::Integer) where {T<:AbstractFloat}
    if n == 0
        return T[]
    elseif n == 1
        return T[one(T)]
    end

    w = ones(T, n)

    @inbounds begin
        half = one(T) / 2
        w[(end - 1):-2:1] .= -1
        w[1] *= half
        w[end] = half
    end

    return w
end

function cheb2_barywts(n::Integer)
    return cheb2_barywts(Float64, n)
end

export cheb2_barywts
