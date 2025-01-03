"""
    cheb2_angles([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}

Compute angles for Chebyshev points of the 2nd kind:
``\\theta_k = \\frac{k\\pi}{n-1}, \\quad k = n-1,\\ldots,0``

# Arguments
- `TR`: Type parameter for the angles (e.g., Float64)
- `n`: Number of points
"""
function cheb2_angles(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    @argcheck n >= 0 "n must be nonnegative"

    if n == 0
        return TR[]
    elseif n == 1
        return TR[convert(TR, π) / 2]
    end

    θ = Array{TR}(undef, n)
    nm1 = n - 1
    pi_over_nm1 = convert(TR, π) / nm1

    @inbounds for k in 0:nm1
        θ[n - k] = k * pi_over_nm1
    end

    return θ
end

function cheb2_angles(n::TI) where {TI<:Integer}
    return cheb2_angles(Float64, n)
end

export cheb2_angles
