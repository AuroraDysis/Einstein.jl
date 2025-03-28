using TestItems

@testitem "cheb1_coeffs2vals" begin
    using LinearAlgebra, Einstein.ChebyshevSuite, Test

    # Set tolerance
    tol = 100 * eps()

    @testset "Single coefficient" begin
        c = [sqrt(2)]
        v = cheb1_coeffs2vals(c)
        @test c ≈ v
    end

    @testset "Even case" begin
        # Simple data (even case)
        c = collect(6.0:-1:1)
        # Exact values
        vTrue = [
            -3 * sqrt(6) / 2 - 5 / sqrt(2) + 2 * sqrt(3) + 7
            4 - sqrt(2) / 2
            -3 * sqrt(6) / 2 + 5 / sqrt(2) - 2 * sqrt(3) + 7
            3 * sqrt(6) / 2 - 5 / sqrt(2) - 2 * sqrt(3) + 7
            4 + sqrt(2) / 2
            3 * sqrt(6) / 2 + 5 / sqrt(2) + 2 * sqrt(3) + 7
        ]

        # Test real branch
        v = cheb1_coeffs2vals(c)
        @test norm(v - vTrue, Inf) < tol
        @test all(iszero, imag.(v))
    end

    @testset "Odd case" begin
        # Simple data (odd case)
        c = collect(5.0:-1:1)
        # Exact values
        vTrue = [
            11 / 2 + sqrt(5) - 2 * sqrt((5 + sqrt(5)) / 2) - sqrt((5 - sqrt(5)) / 2)
            11 / 2 - sqrt(5) - 2 * sqrt((5 - sqrt(5)) / 2) + sqrt((5 + sqrt(5)) / 2)
            3
            11 / 2 - sqrt(5) + 2 * sqrt((5 - sqrt(5)) / 2) - sqrt((5 + sqrt(5)) / 2)
            11 / 2 + sqrt(5) + 2 * sqrt((5 + sqrt(5)) / 2) + sqrt((5 - sqrt(5)) / 2)
        ]

        # Test real branch
        v = cheb1_coeffs2vals(c)
        @test norm(v - vTrue, Inf) < tol
        @test all(iszero, imag.(v))
    end

    @testset "Symmetry preservation" begin
        c = kron(ones(10), Matrix{Float64}(I, 2, 2))
        v1 = cheb1_coeffs2vals(c[:, 1])
        v2 = cheb1_coeffs2vals(c[:, 2])
        @test norm(v1 - reverse(v1), Inf) ≈ 0
        @test norm(v2 + reverse(v2), Inf) ≈ 0
    end

    @testset "Operator style" begin
        n = 100
        coeffs = rand(n)
        op = ChebyshevFirstKindSynthesis{Float64}(n)

        # Test operator call
        vals1 = op(coeffs)
        vals2 = cheb1_coeffs2vals(coeffs)
        @test maximum(abs.(vals1 .- vals2)) < tol

        # Test multiple calls
        for _ in 1:10
            coeffs = rand(n)
            vals1 = op(coeffs)
            vals2 = cheb1_coeffs2vals(coeffs)
            @test maximum(abs.(vals1 .- vals2)) < tol
        end
    end
end
