@testitem "cheb_gauss_vals2coeffs" begin
    using LinearAlgebra

    tol = typetol(Float64)

    function vals2coeffs_via_matrix(v::AbstractVector{T}) where {T<:Number}
        n = length(v)
        A = cheb_gauss_vals2coeffs_matrix(T, n)
        return A * v
    end

    @testset "Even case tests" begin
        v = collect(Float64, 1:6)
        cTrue = [
            7 / 2,
            sqrt(6) / 2 + 5 * sqrt(2) / 6,
            0,
            sqrt(2) / 6,
            0,
            sqrt(6) / 2 - 5 * sqrt(2) / 6,
        ]
        c1 = cheb_gauss_vals2coeffs(v)
        c2 = vals2coeffs_via_matrix(v)
        @test isapprox(c1, cTrue, atol=tol)
        @test isapprox(c2, cTrue, atol=tol)
    end

    @testset "Odd case tests" begin
        v = collect(Float64, 1:5)
        cTrue = [
            3,
            (2 / 5) * (sqrt((5 - sqrt(5)) / 2) + 2 * sqrt((5 + sqrt(5)) / 2)),
            0,
            (2 / 5) * (2 * sqrt((5 - sqrt(5)) / 2) - sqrt((5 + sqrt(5)) / 2)),
            0,
        ]
        c1 = cheb_gauss_vals2coeffs(v)
        c2 = vals2coeffs_via_matrix(v)
        @test isapprox(c1, cTrue, atol=tol)
        @test isapprox(c2, cTrue, atol=tol)
    end

    @testset "Operator style" begin
        n = 100
        vals = rand(n)
        ctx = cheb_gauss_vals2coeffs_context(Float64, n)

        # Test operator call
        coeffs1 = cheb_gauss_vals2coeffs!(ctx, vals)
        coeffs2 = cheb_gauss_vals2coeffs(vals)
        @test isapprox(coeffs1, coeffs2, atol=tol)

        # Test multiple calls
        for _ in 1:10
            vals = rand(n)
            coeffs1 = cheb_gauss_vals2coeffs!(ctx, vals)
            coeffs2 = cheb_gauss_vals2coeffs(vals)
            @test isapprox(coeffs1, coeffs2, atol=tol)
        end
    end
end