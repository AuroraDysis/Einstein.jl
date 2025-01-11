using TestItems

@testitem "cheb2_quadwts" begin
    using LinearAlgebra, Einstein.ChebSuite, Test

    @testset "coefficients" begin
        # Test n=0 case
        @test cheb2_quadwts(0) == Float64[]

        # Test n=1 case
        @test cheb2_quadwts(1) ≈ [2.0]

        # Test n=5 case
        w5 = cheb2_quadwts(5)
        @test w5 ≈ [
            0.0666666666666667,
            0.533333333333333,
            0.800000000000000,
            0.533333333333333,
            0.0666666666666667,
        ]

        w6 = cheb2_quadwts(6)
        @test w6 ≈ [
            0.0400000000000000,
            0.360743041200011,
            0.599256958799989,
            0.599256958799989,
            0.360743041200011,
            0.0400000000000000,
        ]
    end

    @testset "functions" begin
        tol = 100 * eps(Float64)
        n = 32
        x1 = cheb2_pts(n)
        w1 = cheb2_quadwts(n)

        f1 = @. sin(2π * x1)
        If1 = dot(f1, w1)

        @test isapprox(If1, 0.0; atol=tol)

        x2 = cheb2_pts(n, 0.0, 1.0)
        w2 = cheb2_quadwts(n, 0.0, 1.0)

        f2 = @. sin(2π * x2)
        If2 = dot(f2, w2)

        @test isapprox(If2, 0.0; atol=tol)
    end
end
