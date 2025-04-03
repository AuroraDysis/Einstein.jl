@testitem "cheb_series_integrate" begin
    using LinearAlgebra

    f_coeffs = [0.5, 0.0, 0.5]
    If_coeffs = cheb_series_integrate(f_coeffs)
    @test If_coeffs ≈ [1/3, 1/4, 0.0, 1/12]
end
