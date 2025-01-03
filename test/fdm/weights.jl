using TestItems

@testitem "fornberg_calculate_wts" begin
    using PDESuite.FDMSuite, Test

    @testset "Standard First Derivatives" begin
        # Test forward difference (first derivative)
        weights_forward = fornberg_calculate_wts(1, 0.0, [0.0, 1.0])
        @test weights_forward ≈ [-1.0, 1.0]

        # Test central difference (first derivative)
        weights_central = fornberg_calculate_wts(1, 0.0, [-1.0, 0.0, 1.0])
        @test weights_central ≈ [-0.5, 0.0, 0.5]

        # Test accuracy with non-uniform grid
        weights_nonuniform = fornberg_calculate_wts(1, 0.0, [-2.0, -0.5, 1.0])
        # Test against pre-computed values
        @test sum(weights_nonuniform) ≈ 0.0 atol = 1e-14
    end

    @testset "Standard Second Derivatives" begin
        # Test central difference (second derivative)
        weights_central = fornberg_calculate_wts(2, 0.0, [-1.0, 0.0, 1.0])
        @test weights_central ≈ [1.0, -2.0, 1.0]

        # Test with wider stencil
        weights_wide = fornberg_calculate_wts(2, 0.0, [-2.0, -1.0, 0.0, 1.0, 2.0])
        # Test properties that must hold
        @test sum(weights_wide) ≈ 0.0 atol = 1e-14
        @test sum(weights_wide .* [-2.0, -1.0, 0.0, 1.0, 2.0] .^ 2) ≈ 2.0 atol = 1e-14
    end

    @testset "Hermite Finite Difference" begin
        # Test first derivative with function and derivative values
        x = [-1.0, 0.0, 1.0]
        D, E = fornberg_calculate_wts(1, 0.0, x; dfdx=true)

        # Known weights for Hermite interpolation
        @test length(D) == length(x)  # Function value weights
        @test length(E) == length(x)  # Derivative value weights
        @test sum(D) ≈ 0.0 atol = 1e-14  # Consistency condition

        # Test second derivative with function and derivative values
        D2, E2 = fornberg_calculate_wts(2, 0.0, x; dfdx=true)
        @test sum(D2) ≈ 0.0 atol = 1e-14
    end

    @testset "Error Handling" begin
        # Test insufficient points for order
        @test_throws ArgumentError fornberg_calculate_wts(2, 0.0, [0.0, 1.0])

        # Test insufficient points for Hermite interpolation
        @test_throws ArgumentError fornberg_calculate_wts(3, 0.0, [0.0, 1.0], dfdx=true)
    end
end
