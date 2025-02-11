using TestItems

@testitem "cheb_coeffs_eval" begin
    using LinearAlgebra, Einstein.ChebyshevSuite, Test

    tol = 10 * eps()

    @testset "Single coefficient tests" begin
        # Scalar evaluation
        c = [sqrt(2)]
        v = cheb_coeffs_eval(c, 0.0)
        @test v ≈ sqrt(2)

        # Vector evaluation
        x = [-0.5, 1.0]
        v = map(xi -> cheb_coeffs_eval(c, xi), x)
        @test all(v .≈ sqrt(2))
    end

    @testset "Vector coefficient tests" begin
        # Simple polynomial test
        c = collect(5.0:-1:1)
        x = [-0.5, -0.1, 1.0]

        v = map(xi -> cheb_coeffs_eval(c, xi), x)
        v_true = [3.0, 3.1728, 15.0]
        @test norm(v - v_true, Inf) < tol

        # Multiple polynomials
        c2 = reverse(c)
        v = map(xi -> [cheb_coeffs_eval(c, xi), cheb_coeffs_eval(c2, xi)], x)
        v_true = [
            3.0 0.0
            3.1728 3.6480
            15.0 15.0
        ]
        @test norm(hcat(v...)' - v_true, Inf) < tol
    end

    @testset "Edge cases" begin
        # Empty coefficient vector
        @test_throws ArgumentError cheb_coeffs_eval(Float64[], 0.0)

        # Single coefficient
        @test cheb_coeffs_eval([1.0], 0.5) ≈ 1.0
    end
end
