@testitem "cheb_gauss_quadrature_weights" begin
    using LinearAlgebra

    @testset "coefficients" begin
        # Test n=0 case
        @test cheb_gauss_quadrature_weights(0) == Float64[]

        # Test n=1 case
        @test cheb_gauss_quadrature_weights(1) ≈ [2.0]

        # Test n=5 case
        w5 = cheb_gauss_quadrature_weights(5)
        @test w5 ≈ [
            0.167781228466683,
            0.525552104866650,
            0.613333333333333,
            0.525552104866650,
            0.167781228466684,
        ]

        w6 = cheb_gauss_quadrature_weights(6)
        @test w6 ≈ [
            0.118661021381236,
            0.377777777777778,
            0.503561200840986,
            0.503561200840986,
            0.377777777777778,
            0.118661021381236,
        ]
    end

    @testset "functions" begin
        tol = 100 * eps(Float64)
        n = 32
        x1 = cheb_gauss_points(n)
        w1 = cheb_gauss_quadrature_weights(n)

        f1 = @. sin(2π * x1)
        If1 = dot(f1, w1)

        @test isapprox(If1, 0.0; atol=tol)

        x2 = cheb_gauss_points(n, 0.0, 1.0)
        w2 = cheb_gauss_quadrature_weights(n, 1.0)

        f2 = @. sin(2π * x2)
        If2 = dot(f2, w2)

        @test isapprox(If2, 0.0; atol=tol)
    end
end
