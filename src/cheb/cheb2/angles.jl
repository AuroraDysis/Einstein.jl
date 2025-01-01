"""
    cheb2_angles([TR=Float64], n::Integer)

Compute angles for Chebyshev points of the second kind.

# Arguments
- `TR`: Type parameter for the angles (e.g., Float64)
- `n`: Number of points

# Returns
- Vector of n angles in [0,π], ordered decreasing
- Empty vector if n=0
- [π/2] if n=1

# Mathematical Details
For n > 1:
``\\theta_k = \\frac{k\\pi}{n-1}, \\quad k = n-1,\\ldots,0``

These angles generate second-kind Chebyshev points via: x_k = -cos(θ_k)
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

@testset "cheb2_angles" begin
    @testset "n = 5" begin
        n = 5
        θ = cheb2_angles(n)

        @test length(θ) == n
        @test θ ≈
            [3.14159265358979, 2.35619449019235, 1.57079632679490, 0.785398163397448, 0.0]
    end

    @testset "n = 6" begin
        n = 6
        θ = cheb2_angles(n)

        @test length(θ) == n
        @test θ ≈ [
            3.14159265358979,
            2.51327412287183,
            1.88495559215388,
            1.25663706143592,
            0.628318530717959,
            0.0,
        ]
    end
end
