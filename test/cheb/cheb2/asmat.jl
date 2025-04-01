using TestItems

@testitem "vals2coeffs_matrix, coeffs2vals_matrix" begin
    using Einstein.ChebyshevSuite, Test

    n = 32  # Enough points for good accuracy
    x = points(Float64, n)
    A = vals2coeffs_matrix(Float64, n)
    S = coeffs2vals_matrix(Float64, n)

    @testset "Transform and recover" begin
        # Test with polynomial that should be exactly represented
        f = @. 3x^2 + 2x - 1
        coeffs = A * f
        coeffs2 = vals2coeffs(f)
        @test coeffs ≈ coeffs2 rtol = 1e-12
        f_recovered = S * coeffs
        @test f_recovered ≈ f rtol = 1e-12

        # Test with trigonometric function
        f = @. sin(π * x) * cos(2π * x)
        coeffs = A * f
        coeffs2 = vals2coeffs(f)
        @test coeffs ≈ coeffs2 rtol = 1e-12
        f_recovered = S * coeffs
        @test f_recovered ≈ f rtol = 1e-12
    end
end
