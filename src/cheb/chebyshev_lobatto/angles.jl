@doc raw"""
    cheb_lobatto_angles(TF, n) where {TF<:AbstractFloat}

Compute angles for Chebyshev points of the 2nd kind:
```math
\theta_k = \frac{k\pi}{n-1}, \quad k = n-1,\ldots,0
```

# Arguments
- `TF`: Type parameter for the angles (e.g., Float64)
- `n`: Number of points
"""
function cheb_lobatto_angles(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    @argcheck n >= 0 "n must be nonnegative"

    if n == 0
        return TF[]
    elseif n == 1
        return TF[convert(TF, π) / 2]
    end

    θ = Vector{TF}(undef, n)

    nm1 = n - 1
    pi_over_nm1 = convert(TF, π) / nm1
    @inbounds for k in 0:nm1
        θ[n - k] = k * pi_over_nm1
    end

    return θ
end

function cheb_lobatto_angles(n::Integer)
    return cheb_lobatto_angles(Float64, n)
end

export cheb_lobatto_angles
