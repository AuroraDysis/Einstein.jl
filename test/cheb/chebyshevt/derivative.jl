using TestItems

@testitem "chebyshevt_derivative!" begin
    using Einstein.ChebyshevSuite, Test

    @testset "Basic functionality" begin
        # Test case 1: Simple polynomial
        c = [1.0, 2.0, 3.0]
        der = rand(2)
        chebyshevt_derivative!(c, der)
        @test der ≈ [2.0, 12.0]

        # Test case 2: Higher degree
        c = [1.0, 2.0, 3.0, 4.0, 5.0]
        der = rand(6)
        chebyshevt_derivative!(c, der)
        @test der[1:4] ≈ [14.0, 52.0, 24.0, 40.0]
        @test der[5:end] ≈ zeros(2)
    end

    @testset "Edge cases" begin
        # Single coefficient
        c = [1.0]
        der = zeros(0)
        chebyshevt_derivative!(c, der)
        @test isempty(der)

        # Two coefficients
        c = [1.0, 2.0]
        der = zeros(1)
        chebyshevt_derivative!(c, der)
        @test der ≈ [2.0]
    end

    @testset "Known derivatives" begin
        # Test T₃(x) = 4x³ - 3x
        c = [0.0, 3.0, 0.0, 4.0]  # Coefficients of T₃
        der = zeros(3)
        chebyshevt_derivative!(c, der)
        @test der ≈ [15.0, 0.0, 24.0]

        # Test T₄(x) = 8x⁴ - 8x² + 1
        c = [1.0, 0.0, -8.0, 0.0, 8.0]  # Coefficients of T₄
        der = zeros(4)
        chebyshevt_derivative!(c, der)
        @test der ≈ [0.0, 32.0, 0.0, 64.0]
    end
end
