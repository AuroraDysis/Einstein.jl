@doc raw"""
    cheb2_angles([T=Float64], n::Integer) where {T<:AbstractFloat}

Compute angles for Chebyshev points of the 2nd kind:
```math
\theta_k = \frac{k\pi}{n-1}, \quad k = n-1,\ldots,0
```

# Arguments
- `T`: Type parameter for the angles (e.g., Float64)
- `n`: Number of points
"""
function cheb2_angles(::Type{T}, n::Integer) where {T<:AbstractFloat}
    @argcheck n >= 0 "n must be nonnegative"

    if n == 0
        return T[]
    elseif n == 1
        return T[convert(T, π) / 2]
    end

    θ = Array{T}(undef, n)
    nm1 = n - 1
    pi_over_nm1 = convert(T, π) / nm1

    @inbounds for k in 0:nm1
        θ[n - k] = k * pi_over_nm1
    end

    return θ
end

function cheb2_angles(n::Integer)
    return cheb2_angles(Float64, n)
end

export cheb2_angles
