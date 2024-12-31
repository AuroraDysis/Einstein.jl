using Test
using GRSuite.ChebyshevSuite

@testset "Chebyshev Spectral Transforms" begin
    n = 32  # Enough points for good accuracy
    x = cheb2_grid(Float64, n)
    A, S = cheb2_asmat(Float64, n)

    @testset "Transform and recover" begin
        # Test with polynomial that should be exactly represented
        f = @. 3x^2 + 2x - 1
        coeffs = A * f
        f_recovered = S * coeffs
        @test f_recovered ≈ f rtol = 1e-12

        # Test with trigonometric function
        f = @. sin(π * x) * cos(2π * x)
        coeffs = A * f
        f_recovered = S * coeffs
        @test f_recovered ≈ f rtol = 1e-12
    end
end

@testset "Chebyshev Integration" begin
    @testset "Standard domain [-1,1]" begin
        n = 32
        x = cheb2_grid(Float64, n)
        intmat = cheb2_intmat(Float64, n)

        # Test 1: Polynomial integration
        f = @. x^3  # f(x) = x³
        F_numeric = intmat * f   # Should give (x⁴ - 1) / 4
        F_exact = @. (x^4 - 1) / 4
        @test F_numeric ≈ F_exact rtol = 1e-12

        # Test 2: Trigonometric integration
        f = @. sin(π * x)
        F_numeric = intmat * f   # Should give -(cos(πx)+1)/π
        F_exact = @. -(cos(π * x) + 1) / π
        @test F_numeric ≈ F_exact rtol = 1e-12
    end

    @testset "Mapped domain [0,π]" begin
        n = 32
        intmat = cheb2_intmat(Float64, n, 0.0, Float64(π))
        x = cheb2_grid(Float64, n, 0.0, Float64(π))

        # Test: Integration of sin(x) from 0 to x
        f = sin.(x)
        F_numeric = intmat * f   # Should give -cos(x) + 1
        F_exact = @. -cos(x) + 1
        @test F_numeric ≈ F_exact rtol = 1e-12

        # Test: Integration of x*cos(x)
        f = @. x * cos(x)
        F_numeric = intmat * f   # Should give x*sin(x) - sin(x)
        F_exact = @. x * sin(x) + cos(x) - 1
        @test F_numeric ≈ F_exact rtol = 1e-12
    end
end
