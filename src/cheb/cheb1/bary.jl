"""
    cheb1_bary_wts([TR=Float64], n::Integer)

Compute the barycentric weights for Chebyshev points of the first kind.

# Arguments
- `TR`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# Returns
- Vector of n barycentric weights

# Mathematical Details
For Chebyshev points of the first kind, the barycentric weights are:
```math
w_j = (-1)^j \\sin\\left(\\frac{(2j+1)\\pi}{2n}\\right), \\quad j = 0,\\ldots,n-1
```

These weights provide optimal stability for barycentric interpolation at
Chebyshev points of the first kind.

See also: [`bary`](@ref), [`cheb1_grid`](@ref)
"""
function cheb1_bary_wts(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    if n == 0
        return TR[]
    elseif n == 1
        return [one(TR)]
    end

    w = zeros(TR, n)
    pi_over_n = convert(TR, π) / n
    
    @inbounds for j in 0:(n-1)
        θ = (n - j - 0.5) * pi_over_n
        w[j + 1] = sin(θ)
    end

    # Flip signs for alternate elements
    @inbounds @. w[(end-1):-2:1] *= -1

    return w
end

function cheb1_bary_wts(n::TI) where {TI<:Integer}
    return cheb1_bary_wts(Float64, n)
end

export cheb1_bary_wts
