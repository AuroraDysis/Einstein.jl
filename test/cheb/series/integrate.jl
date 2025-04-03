@testitem "cheb_series_integrate" begin
    using LinearAlgebra

    @testset "Simple polynomial f'(x) = 1 + x^2" begin
        f_coeffs = [0.5, 0.0, 0.5]
        If_coeffs = cheb_series_integrate(f_coeffs)
        @test If_coeffs ≈ [1/3, 1/4, 0.0, 1/12]
    end

    @testset "Constant function f'(x) = 1" begin
        f_coeffs = [1.0]
        If_coeffs = cheb_series_integrate(f_coeffs)
        @test If_coeffs ≈ [1.0, 1.0]
    end

    @testset "Linear function f'(x) = x" begin
        f_coeffs = [0.0, 1.0]
        If_coeffs = cheb_series_integrate(f_coeffs)
        @test If_coeffs ≈ [-0.25, 0.0, 0.25]
    end

    @testset "Higher degree polynomial f'(x) = x^3" begin
        f_coeffs = [0.0, 0.0, 0.0, 1.0]
        If_coeffs = cheb_series_integrate(f_coeffs)
        @test If_coeffs ≈ [0.125, 0.0, -0.25, 0.0, 0.125]
    end

    @testset "In-place integration" begin
        f_coeffs = [0.5, 0.0, 0.5]
        result = zeros(4)
        cheb_series_integrate!(result, f_coeffs)
        @test result ≈ [1/3, 1/4, 0.0, 1/12]
    end

    @testset "Zero function f'(x) = 0" begin
        f_coeffs = [0.0, 0.0, 0.0]
        If_coeffs = cheb_series_integrate(f_coeffs)
        @test If_coeffs ≈ [0.0, 0.0, 0.0, 0.0]
    end

    @testset "Negative coefficients" begin
        f_coeffs = [-0.5, 0.0, -0.5]
        If_coeffs = cheb_series_integrate(f_coeffs)
        @test If_coeffs ≈ [-1/3, -1/4, 0.0, -1/12]
    end

    @testset "Mixed coefficients" begin
        f_coeffs = [1.0, -0.5, 0.25, -0.125]
        If_coeffs = cheb_series_integrate(f_coeffs)
        @test length(If_coeffs) == 5  # Check output length
        @test If_coeffs[1] ≈ sum((-1)^(r+1) * If_coeffs[r+1] for r in 1:4)  # Check constant of integration
    end

    @testset "Complex coefficients" begin
        f_coeffs = [1.0 + 2.0im, 3.0 - 1.0im, 2.0 + 0.0im]
        If_coeffs = cheb_series_integrate(f_coeffs)
        @test If_coeffs ≈ [1.0 + 2.0im, 3.0 - 1.0im, 2.0 + 0.0im, 0.0, 0.0]
    end
end
