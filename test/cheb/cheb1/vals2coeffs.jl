using TestItems

@testitem "cheb1_vals2coeffs" begin
    using LinearAlgebra, Einstein.ChebyshevSuite, Test

    tol = 100 * eps()

    @testset "Single value conversion" begin
        v = sqrt(2)
        c = cheb1_vals2coeffs([v])
        @test v ≈ c[1]
    end

    @testset "Even case tests" begin
        v = Float64[1:6;]
        cTrue = [
            7 / 2,
            sqrt(6) / 2 + 5 * sqrt(2) / 6,
            0,
            sqrt(2) / 6,
            0,
            sqrt(6) / 2 - 5 * sqrt(2) / 6,
        ]
        c = cheb1_vals2coeffs(v)
        @test norm(c - cTrue, Inf) < tol
        @test all(x -> abs(imag(x)) < tol, c)
    end

    @testset "Odd case tests" begin
        v = Float64[1:5;]
        cTrue = [
            3,
            (2 / 5) * (sqrt((5 - sqrt(5)) / 2) + 2 * sqrt((5 + sqrt(5)) / 2)),
            0,
            (2 / 5) * (2 * sqrt((5 - sqrt(5)) / 2) - sqrt((5 + sqrt(5)) / 2)),
            0,
        ]
        c = cheb1_vals2coeffs(v)
        @test norm(c - cTrue, Inf) < tol
        @test all(x -> abs(imag(x)) < tol, c)
    end

    @testset "Operator style" begin
        n = 100
        vals = rand(n)
        op = Cheb1Vals2CoeffsOp(n)

        # Test operator call
        coeffs1 = op(vals)
        coeffs2 = cheb1_vals2coeffs(vals)
        @test maximum(abs.(coeffs1 .- coeffs2)) < tol

        # Test multiple calls
        for _ in 1:10
            vals = rand(n)
            coeffs1 = op(vals)
            coeffs2 = cheb1_vals2coeffs(vals)
            @test maximum(abs.(coeffs1 .- coeffs2)) < tol
        end
    end
end