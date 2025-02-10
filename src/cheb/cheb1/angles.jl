@doc raw"""
    cheb1_angles([TF=Float64], n::Integer) where {TF<:AbstractFloat}
    cheb1_angles!(θ::Vector{TF}, n::Integer) where {TF<:AbstractFloat}

Compute angles for Chebyshev points of the 1st kind:
```math
\theta_k = \frac{(2k + 1)\pi}{2n}, \quad k = n-1,\ldots,0
```

# Arguments
- `TF`: Type parameter for the angles (e.g., Float64)
- `n`: Number of points
"""
function cheb1_angles(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    @argcheck n >= 0 "n must be nonnegative"

    θ = Array{TF}(undef, n)
    cheb1_angles!(θ, n)

    return θ
end

function cheb1_angles(n::Integer)
    return cheb1_angles(Float64, n)
end

function cheb1_angles!(θ::Vector{TF}, n::Integer) where {TF<:AbstractFloat}
    @argcheck length(θ) == n "length(θ) must be equal to n"

    if n == 1
        θ[1] = convert(TF, π) / 2
        return nothing
    end

    pi_over_2n = convert(TF, pi) / (2 * n)
    @inbounds for k in 0:(n - 1)
        θ[n - k] = (2 * k + 1) * pi_over_2n
    end

    return nothing
end

export cheb1_angles, cheb1_angles!
