@testset "cheb1_amat, cheb1_smat" begin
    n = 32  # Enough points for good accuracy
    x = cheb1_pts(Float64, n)
    A = cheb1_amat(Float64, n)
    S = cheb1_smat(Float64, n)

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
