"""
    cheb1_angles([TR=Float64], n::Integer)

Compute angles for Chebyshev points of the first kind:
``\\theta_k = \\frac{(2k + 1)\\pi}{2n}, \\quad k = n-1,\\ldots,0``

# Arguments
- `TR`: Type parameter for the angles (e.g., Float64)
- `n`: Number of points

# Returns
- Vector of n angles in [0,π], ordered decreasing
"""
function cheb1_angles(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    @argcheck n >= 0 "n must be nonnegative"

    if n == 0
        return TR[]
    elseif n == 1
        return TR[convert(TR, π) / 2]
    end

    θ = Array{TR}(undef, n)

    pi_over_2n = convert(TR, pi) / (2 * n)
    @inbounds for k in 0:(n - 1)
        θ[n - k] = (2 * k + 1) * pi_over_2n
    end

    return θ
end

function cheb1_angles(n::TI) where {TI<:Integer}
    return cheb1_angles(Float64, n)
end

export cheb1_angles

@testset "cheb1_angles" begin
    @testset "n = 5" begin
        n = 5
        θ = cheb1_angles(n)

        @test length(θ) == n
        @test θ ≈ [
            2.82743338823081,
            2.19911485751286,
            1.57079632679490,
            0.942477796076938,
            0.314159265358979,
        ]
    end

    @testset "n = 6" begin
        n = 6
        θ = cheb1_angles(n)

        @test length(θ) == n
        @test θ ≈ [
            2.87979326579064,
            2.35619449019235,
            1.83259571459405,
            1.30899693899575,
            0.785398163397448,
            0.261799387799149,
        ]
    end
end
