@doc raw"""
    angles(TF, n) where {TF<:AbstractFloat}

Compute angles for Chebyshev points of the 1st kind:
```math
\theta_k = \frac{(2k + 1)\pi}{2n}, \quad k = n-1,\ldots,0
```

# Arguments
- `TF`: Type parameter for the angles (e.g., Float64)
- `n`: Number of points
"""
function angles(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    @argcheck n >= 0 "n must be nonnegative"

    if n == 0
        return TF[]
    elseif n == 1
        return TF[convert(TF, π) / 2]
    end

    θ = Vector{TF}(undef, n)
    pi_over_2n = convert(TF, pi) / (2 * n)
    @inbounds for k in 0:(n - 1)
        θ[n - k] = (2 * k + 1) * pi_over_2n
    end

    return θ
end

function angles(n::Integer)
    return angles(Float64, n)
end

export angles
