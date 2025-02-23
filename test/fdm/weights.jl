using TestItems

@testitem "fdm_weights_fornberg" begin
    using Einstein.FiniteDifferenceSuite, Test

    @testset "Standard First Derivatives" begin
        # Test forward difference (first derivative)
        weights_forward = fdm_weights_fornberg(1, 0.0, [0.0, 1.0])
        @test weights_forward ≈ [-1.0, 1.0]

        # Test central difference (first derivative)
        weights_central = fdm_weights_fornberg(1, 0.0, [-1.0, 0.0, 1.0])
        @test weights_central ≈ [-0.5, 0.0, 0.5]

        # Test accuracy with non-uniform grid
        weights_nonuniform = fdm_weights_fornberg(1, 0.0, [-2.0, -0.5, 1.0])
        # Test against pre-computed values
        @test sum(weights_nonuniform) ≈ 0.0 atol = 1e-14
    end

    @testset "Standard Second Derivatives" begin
        # Test central difference (second derivative)
        weights_central = fdm_weights_fornberg(2, 0.0, [-1.0, 0.0, 1.0])
        @test weights_central ≈ [1.0, -2.0, 1.0]

        # Test with wider stencil
        weights_wide = fdm_weights_fornberg(2, 0.0, [-2.0, -1.0, 0.0, 1.0, 2.0])
        # Test properties that must hold
        @test sum(weights_wide) ≈ 0.0 atol = 1e-14
        @test sum(weights_wide .* [-2.0, -1.0, 0.0, 1.0, 2.0] .^ 2) ≈ 2.0 atol = 1e-14
    end

    @testset "Hermite Finite Difference" begin
        # Test first derivative with function and derivative values
        x = [-1.0, 0.0, 1.0]
        D, E = fdm_weights_fornberg(1, 0.0, x; hermite=true)

        # Known weights for Hermite interpolation
        @test length(D) == length(x)  # Function value weights
        @test length(E) == length(x)  # Derivative value weights
        @test sum(D) ≈ 0.0 atol = 1e-14  # Consistency condition

        # Test second derivative with function and derivative values
        D2, E2 = fdm_weights_fornberg(2, 0.0, x; hermite=true)
        @test sum(D2) ≈ 0.0 atol = 1e-14
    end

    @testset "Error Handling" begin
        # Test insufficient points for order
        @test_throws ArgumentError fdm_weights_fornberg(2, 0.0, [0.0, 1.0])

        # Test insufficient points for Hermite interpolation
        @test_throws ArgumentError fdm_weights_fornberg(3, 0.0, [0.0, 1.0], hermite=true)
    end
end

@testitem "fdm_central_weights" begin
    using Einstein.FiniteDifferenceSuite, Test

    # Test 1st derivative, 2nd order accuracy
    stencil = [-1//2, 0//1, 1//2]
    @test fdm_central_weights(1, 2) == stencil
    for type in [Float64, BigFloat]
        @test fdm_central_weights(type, 1, 2) ≈ type.(stencil)
    end

    # Test 2nd derivative, 2nd order accuracy 
    stencil = [1//1, -2//1, 1//1]
    @test fdm_central_weights(2, 2) ≈ stencil
    for type in [Float64, BigFloat]
        @test fdm_central_weights(type, 2, 2) ≈ type.(stencil)
    end

    # Test 1st derivative, 4th order accuracy
    stencil = [1//12, -2//3, 0//1, 2//3, -1//12]
    @test fdm_central_weights(1, 4) == stencil
    for type in [Float64, BigFloat]
        @test fdm_central_weights(type, 1, 4) ≈ type.(stencil)
    end

    # Test 6th derivative, 6th order accuracy
    stencil = [
        13//240,
        -19//24,
        87//16,
        -39//2,
        323//8,
        -1023//20,
        323//8,
        -39//2,
        87//16,
        -19//24,
        13//240,
    ]
    @test fdm_central_weights(6, 6) ≈ stencil
    for type in [Float64, BigFloat]
        @test fdm_central_weights(type, 6, 6) ≈ type.(stencil)
    end
end

@testitem "fdm_hermite_weights" begin
    using Einstein.FiniteDifferenceSuite, Test

    # Test 2nd derivative, 4th order accuracy
    D2, E2 = fdm_hermite_weights(2, 4)
    @test D2 == [2//1, -4//1, 2//1]
    @test E2 == [1//2, 0//1, -1//2]

    # Test 2th derivative, 8nd order accuracy
    D4, E4 = fdm_hermite_weights(2, 8)
    @test D4 == [7//54, 64//27, -5//1, 64//27, 7//54]
    @test E4 == [1//36, 8//9, 0//1, -8//9, -1//36]
end

@testitem "fdm_extrapwts" begin
    using Einstein.FiniteDifferenceSuite, Test

    @test fdm_extrapwts_left(4) == [4, -6, 4, -1]
    @test fdm_extrapwts_right(4) == [-1, 4, -6, 4]

    @test fdm_extrapwts_left(5) == [5, -10, 10, -5, 1]
    @test fdm_extrapwts_right(5) == [1, -5, 10, -10, 5]
end

@testitem "fdm_boundary_weights" begin
    using Einstein.FiniteDifferenceSuite, Test

    Dl14, Dr14 = fdm_boundary_weights(1, 4)
    @test Dl14[:, 1] == [-25//12, 4//1, -3//1, 4//3, -1//4]
    @test Dl14[:, 2] == [-1//4, -5//6, 3//2, -1//2, 1//12]
    @test Dr14[:, 1] == [-1//12, 1//2, -3//2, 5//6, 1//4]
    @test Dr14[:, 2] == [1//4, -4//3, 3//1, -4//1, 25//12]

    Dl24, Dr24 = fdm_boundary_weights(2, 4)
    @test Dl24[:, 1] == [15//4, -77//6, 107//6, -13//1, 61//12, -5//6]
    @test Dl24[:, 2] == [5//6, -5//4, -1//3, 7//6, -1//2, 1//12]
    @test Dr24[:, 1] == [1//12, -1//2, 7//6, -1//3, -5//4, 5//6]
    @test Dr24[:, 2] == [-5//6, 61//12, -13//1, 107//6, -77//6, 15//4]
end
