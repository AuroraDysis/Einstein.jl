using TestItems

@testitem "cheb2_analysis_matrix, cheb2_synthesis_matrix" begin
    using Einstein.ChebyshevSuite, Test

    n = 32  # Enough points for good accuracy
    x = cheb2_points(Float64, n)
    A = cheb2_analysis_matrix(Float64, n)
    S = cheb2_synthesis_matrix(Float64, n)

    @testset "Transform and recover" begin
        # Test with polynomial that should be exactly represented
        f = @. 3x^2 + 2x - 1
        coeffs = A * f
        coeffs2 = cheb2_vals2coeffs(f)
        @test coeffs ≈ coeffs2 rtol = 1e-12
        f_recovered = S * coeffs
        @test f_recovered ≈ f rtol = 1e-12

        # Test with trigonometric function
        f = @. sin(π * x) * cos(2π * x)
        coeffs = A * f
        coeffs2 = cheb2_vals2coeffs(f)
        @test coeffs ≈ coeffs2 rtol = 1e-12
        f_recovered = S * coeffs
        @test f_recovered ≈ f rtol = 1e-12
    end
end
