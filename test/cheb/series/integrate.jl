@testitem "cheb_series_integrate" begin
    using LinearAlgebra

    # Test case 1: Simple polynomial f'(x) = 1 + x^2
    f_coeffs = [0.5, 0.0, 0.5]
    If_coeffs = cheb_series_integrate(f_coeffs)
    @test If_coeffs ≈ [1/3, 1/4, 0.0, 1/12]

    # Test case 2: Constant function f'(x) = 1
    f_coeffs = [1.0]
    If_coeffs = cheb_series_integrate(f_coeffs)
    @test If_coeffs ≈ [1.0, 1.0]

    # Test case 3: Linear function f'(x) = x
    f_coeffs = [0.0, 1.0]
    If_coeffs = cheb_series_integrate(f_coeffs)
    @test If_coeffs ≈ [-0.25, 0.0, 0.25]

    # Test case 4: Higher degree polynomial f'(x) = x^3
    f_coeffs = [0.0, 0.0, 0.0, 1.0]
    If_coeffs = cheb_series_integrate(f_coeffs)
    @test If_coeffs ≈ [0.125, 0.0, -0.25, 0.0, 0.125]

    # Test case 5: In-place integration
    f_coeffs = [0.5, 0.0, 0.5]
    result = zeros(4)
    cheb_series_integrate!(result, f_coeffs)
    @test result ≈ [1/3, 1/4, 0.0, 1/12]

    # Test case 6: Zero function f'(x) = 0
    f_coeffs = [0.0, 0.0, 0.0]
    If_coeffs = cheb_series_integrate(f_coeffs)
    @test If_coeffs ≈ [0.0, 0.0, 0.0, 0.0]

    # Test case 7: Negative coefficients
    f_coeffs = [-0.5, 0.0, -0.5]
    If_coeffs = cheb_series_integrate(f_coeffs)
    @test If_coeffs ≈ [-1/3, -1/4, 0.0, -1/12]

    # Test case 8: Mixed coefficients
    f_coeffs = [1.0, -0.5, 0.25, -0.125]
    If_coeffs = cheb_series_integrate(f_coeffs)
    @test length(If_coeffs) == 5  # Check output length
    @test If_coeffs[1] ≈ sum((-1)^(r+1) * If_coeffs[r+1] for r in 1:4)  # Check constant of integration
end
