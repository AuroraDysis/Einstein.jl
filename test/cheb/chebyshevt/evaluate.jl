using TestItems

@testitem "chebyshevt_evaluate" begin
    using LinearAlgebra, Einstein.ChebyshevSuite, Test

    tol = 10 * eps()

    @testset "Single coefficient tests" begin
        # Scalar evaluation
        c = [sqrt(2)]
        v = chebyshevt_evaluate(c, 0.0)
        @test v ≈ sqrt(2)

        # Vector evaluation
        x = [-0.5, 1.0]
        v = chebyshevt_evaluate(c, x)
        @test all(v .≈ sqrt(2))
    end

    @testset "Vector coefficient tests" begin
        # Simple polynomial test
        c = collect(5.0:-1:1)
        x = [-0.5, -0.1, 1.0]

        v = chebyshevt_evaluate(c, x)
        v_true = [3.0, 3.1728, 15.0]
        @test norm(v - v_true, Inf) < tol

        # Multiple polynomials
        c2 = reverse(c)
        v = map(xi -> [chebyshevt_evaluate(c, xi), chebyshevt_evaluate(c2, xi)], x)
        v_true = [
            3.0 0.0
            3.1728 3.6480
            15.0 15.0
        ]
        @test norm(hcat(v...)' - v_true, Inf) < tol
    end

    @testset "Complex coefficients" begin
        # Test with complex coefficients
        c = [1.0 + 2.0im, 3.0 - 1.0im, 2.0 + 0.0im]
        x = [0.5, -0.5]
        
        v = chebyshevt_evaluate(c, x)
        # Expected values calculated using direct evaluation:
        # T₀(x) = 1
        # T₁(x) = x
        # T₂(x) = 2x² - 1
        v_true = [
            (1.0 + 2.0im) + (3.0 - 1.0im) * 0.5 + (2.0 + 0.0im) * (2 * 0.5^2 - 1),
            (1.0 + 2.0im) + (3.0 - 1.0im) * (-0.5) + (2.0 + 0.0im) * (2 * (-0.5)^2 - 1)
        ]
        @test norm(v - v_true, Inf) < tol
    end

    @testset "Edge cases" begin
        # Empty coefficient vector
        @test_throws ArgumentError chebyshevt_evaluate(Float64[], 0.0)

        # Single coefficient
        @test chebyshevt_evaluate([1.0], 0.5) ≈ 1.0
    end
end
