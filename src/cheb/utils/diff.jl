"""
    cheb_diff(coeffs::VT) where {T<:AbstractFloat,VT<:AbstractVector{T}}
    cheb_diff!(coeffs::VT1, coeffs_der::VT2) where {T<:AbstractFloat,VT1<:AbstractVector{T},VT2<:AbstractVector{T}}

Compute derivatives of Chebyshev coefficients.

# Arguments
- `coeffs`: Input vector of Chebyshev coefficients with length n
- `coeffs_der`: Pre-allocated output vector for derivative coefficients (length at least n - 1)
"""
function cheb_diff!(
    coeffs::VT1, coeffs_der::VT2
) where {T<:AbstractFloat,VT1<:AbstractVector{T},VT2<:AbstractVector{T}}
    n = length(coeffs)
    n_der = length(coeffs_der)

    @argcheck n >= 1 "coeffs must have at least one element"
    @argcheck n_der >= n - 1 "coeffs_der must have at least n - 1 elements"

    if n == 1
        coeffs_der .= 0
        return nothing
    end

    # Compute 2k * c[k+1] for k = 1:n-1
    for i in 1:(n - 1)
        coeffs_der[i] = 2 * i * coeffs[i + 1]
    end
    coeffs_der[n:end] .= 0

    # Compute cumulative sums for odd and even indices separately
    # This avoids the need to create temporary arrays
    begin
        # Odd indices (n-1:-2:1)
        for i in (n - 1):-2:3
            coeffs_der[i - 2] += coeffs_der[i]
        end

        # Even indices (n-2:-2:1)
        for i in (n - 2):-2:3
            coeffs_der[i - 2] += coeffs_der[i]
        end

        # Scale first coefficient
        half = one(T) / 2
        coeffs_der[1] *= half
    end

    return nothing
end

function cheb_diff(coeffs::VT) where {T<:AbstractFloat,VT<:AbstractVector{T}}
    n = length(coeffs)
    coeffs_der = similar(coeffs, n - 1)
    cheb_diff!(coeffs, coeffs_der)
    return coeffs_der
end

@testset "cheb_diff!" begin
    @testset "Basic functionality" begin
        # Test case 1: Simple polynomial
        c = [1.0, 2.0, 3.0]
        der = rand(2)
        cheb_diff!(c, der)
        @test der ≈ [2.0, 12.0]

        # Test case 2: Higher degree
        c = [1.0, 2.0, 3.0, 4.0, 5.0]
        der = rand(6)
        cheb_diff!(c, der)
        @test der[1:4] ≈ [14.0, 52.0, 24.0, 40.0]
        @test der[5:end] ≈ zeros(2)
    end

    @testset "Edge cases" begin
        # Single coefficient
        c = [1.0]
        der = zeros(0)
        cheb_diff!(c, der)
        @test isempty(der)

        # Two coefficients
        c = [1.0, 2.0]
        der = zeros(1)
        cheb_diff!(c, der)
        @test der ≈ [2.0]
    end

    @testset "Known derivatives" begin
        # Test T₃(x) = 4x³ - 3x
        c = [0.0, 3.0, 0.0, 4.0]  # Coefficients of T₃
        der = zeros(3)
        cheb_diff!(c, der)
        @test der ≈ [15.0, 0.0, 24.0]

        # Test T₄(x) = 8x⁴ - 8x² + 1
        c = [1.0, 0.0, -8.0, 0.0, 8.0]  # Coefficients of T₄
        der = zeros(4)
        cheb_diff!(c, der)
        @test der ≈ [0.0, 32.0, 0.0, 64.0]
    end
end

export cheb_diff, cheb_diff!
