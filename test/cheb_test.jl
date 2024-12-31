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

@testset "Barycentric Interpolation" begin
    @testset "Barycentric weights" begin
        # Test n=1 case
        w1 = cheb2_interp_wts(Float64, 1)
        @test w1 ≈ [1.0]

        # Test n=2 case
        w2 = cheb2_interp_wts(Float64, 2)
        @test w2 ≈ [-0.5, 0.5]

        # Test n=5 case
        w5 = cheb2_interp_wts(Float64, 5)
        @test w5 ≈ [0.5, -1.0, 1.0, -1.0, 0.5]
    end

    @testset "Polynomial interpolation" begin
        n = 5
        x = cheb2_grid(Float64, n)
        w = cheb2_interp_wts(Float64, n)

        # Test exact interpolation of quadratic polynomial
        f = @. 2x^2 - x + 1
        x0 = 0.3  # Test point
        exact = 2x0^2 - x0 + 1
        interp = cheb2_interp(w, x, f, x0)
        @test interp ≈ exact rtol = 1e-14

        # Test interpolation at nodes (should be exact)
        for i in 1:n
            @test cheb2_interp(w, x, f, x[i]) ≈ f[i] rtol = 1e-14
        end
    end

    @testset "Trigonometric function interpolation" begin
        n = 32  # More points for better accuracy
        x = cheb2_grid(Float64, n)
        w = cheb2_interp_wts(Float64, n)

        # Test interpolation of sin(πx)
        f = @. sin(π * x)
        test_points = [-0.7, -0.3, 0.0, 0.4, 0.8]
        f_exact = sin.(π * test_points)
        f_interp = [cheb2_interp(w, x, f, x0) for x0 in test_points]
        @test f_interp ≈ f_exact rtol = 1e-12
    end
end
