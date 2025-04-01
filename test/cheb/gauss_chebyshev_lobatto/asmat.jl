using TestItems

@testitem "vals2coeffs_matrix, coeffs2vals_matrix" begin
    using Einstein.ChebyshevSuite, Test

    n = 32  # Enough points for good accuracy
    x = GaussChebyshevLobattoGrid.points(Float64, n)
    A = GaussChebyshevLobattoGrid.vals2coeffs_matrix(Float64, n)
    S = GaussChebyshevLobattoGrid.coeffs2vals_matrix(Float64, n)

    tol = typetol(Float64)

    @testset "Transform and recover" begin
        # Test with polynomial that should be exactly represented
        f = @. 3x^2 + 2x - 1
        coeffs = A * f
        coeffs2 = GaussChebyshevLobattoGrid.vals2coeffs(f)
        @test isapprox(coeffs, coeffs2, rtol=tol)
        f_recovered = S * coeffs
        @test isapprox(f_recovered, f, rtol=tol)

        # Test with trigonometric function
        f = @. sin(π * x) * cos(2π * x)
        coeffs = A * f
        coeffs2 = GaussChebyshevLobattoGrid.vals2coeffs(f)
        @test isapprox(coeffs, coeffs2, rtol=tol)
        f_recovered = S * coeffs
        @test isapprox(f_recovered, f, rtol=tol)
    end
end
