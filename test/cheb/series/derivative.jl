@testitem "cheb_series_derivative!" begin
    @testset "In-place version" begin
        # Test case 1: Simple polynomial
        c = [1.0, 2.0, 3.0]
        cheb_series_derivative!(c)
        @test c[1:2] ≈ [2.0, 12.0]
        @test c[3] ≈ 0.0

        # Test case 2: Higher degree
        c = [1.0, 2.0, 3.0, 4.0, 5.0]
        cheb_series_derivative!(c)
        @test c[1:4] ≈ [14.0, 52.0, 24.0, 40.0]
        @test c[5] ≈ 0.0
    end

    @testset "Non-in-place version" begin
        c = [1.0, 2.0, 3.0]
        der = Array{Float64}(undef, 2)
        cheb_series_derivative!(der, c)
        @test der ≈ [2.0, 12.0]

        c = [1.0, 2.0, 3.0, 4.0, 5.0]
        der = Array{Float64}(undef, 4)
        cheb_series_derivative!(der, c)
        @test der ≈ [14.0, 52.0, 24.0, 40.0]

        c = [1.0, 2.0, 3.0]
        der = Array{Float64}(undef, 5)
        cheb_series_derivative!(der, c)
        @test der[1:2] ≈ [2.0, 12.0]
        @test all(isapprox.(der[3:end], 0.0))
    end

    @testset "Copy version" begin
        c = [1.0, 2.0, 3.0]
        der = cheb_series_derivative(c)
        @test der ≈ [2.0, 12.0]

        c = [1.0, 2.0, 3.0, 4.0, 5.0]
        der = cheb_series_derivative(c)
        @test der ≈ [14.0, 52.0, 24.0, 40.0]
    end

    @testset "Edge cases" begin
        # Single coefficient
        c = [1.0]
        cheb_series_derivative!(c)
        @test c[1] ≈ 0.0

        # Two coefficients
        c = [1.0, 2.0]
        cheb_series_derivative!(c)
        @test c[1] ≈ 2.0
        @test c[2] ≈ 0.0
    end

    @testset "Known derivatives" begin
        # Test T₃(x) = 4x³ - 3x
        c = [0.0, 3.0, 0.0, 4.0]  # Coefficients of T₃
        cheb_series_derivative!(c)
        @test c[1:3] ≈ [15.0, 0.0, 24.0]
        @test c[4] ≈ 0.0

        # Test T₄(x) = 8x⁴ - 8x² + 1
        c = [1.0, 0.0, -8.0, 0.0, 8.0]  # Coefficients of T₄
        cheb_series_derivative!(c)
        @test c[1:4] ≈ [0.0, 32.0, 0.0, 64.0]
        @test c[5] ≈ 0.0
    end

    @testset "Complex coefficients" begin
        # Test case 1: Simple polynomial
        c = [1.0im, 2.0im, 3.0im]
        cheb_series_derivative!(c)
        @test c[1:2] ≈ [2.0im, 12.0im]
        @test c[3] ≈ 0.0

        # Test case 2: Higher degree
        c = [1.0im, 2.0im, 3.0im, 4.0im, 5.0im]
        cheb_series_derivative!(c)
        @test c[1:4] ≈ [14.0im, 52.0im, 24.0im, 40.0im]
        @test c[5] ≈ 0.0
    end
end
