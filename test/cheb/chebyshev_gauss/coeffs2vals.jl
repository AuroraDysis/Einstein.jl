@testitem "cheb_gauss_coeffs2vals" begin
    using LinearAlgebra

    # Set tolerance
    tol = typetol(Float64)

    function coeffs2vals_via_matrix(c::AbstractVector{T}) where {T<:Number}
        n = length(c)
        A = cheb_gauss_coeffs2vals_matrix(T, n)
        return A * c
    end

    @testset "Even case" begin
        # Simple data (even case)
        c = collect(Float64, 6:-1:1)
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
        v1 = cheb_gauss_coeffs2vals(c)
        @test isapprox(v1, vTrue, atol=tol)
        @test all(iszero, imag.(v1))

        v2 = coeffs2vals_via_matrix(c)
        @test isapprox(v2, vTrue, atol=tol)
    end

    @testset "Odd case" begin
        # Simple data (odd case)
        c = collect(Float64, 5:-1:1)
        # Exact values
        vTrue = [
            11 / 2 + sqrt(5) - 2 * sqrt((5 + sqrt(5)) / 2) - sqrt((5 - sqrt(5)) / 2)
            11 / 2 - sqrt(5) - 2 * sqrt((5 - sqrt(5)) / 2) + sqrt((5 + sqrt(5)) / 2)
            3
            11 / 2 - sqrt(5) + 2 * sqrt((5 - sqrt(5)) / 2) - sqrt((5 + sqrt(5)) / 2)
            11 / 2 + sqrt(5) + 2 * sqrt((5 + sqrt(5)) / 2) + sqrt((5 - sqrt(5)) / 2)
        ]

        # Test real branch
        v1 = cheb_gauss_coeffs2vals(c)
        @test isapprox(v1, vTrue, atol=tol)
        @test all(iszero, imag.(v1))

        v2 = coeffs2vals_via_matrix(c)
        @test isapprox(v2, vTrue, atol=tol)
    end

    @testset "Symmetry preservation" begin
        c = kron(ones(10), Matrix{Float64}(I, 2, 2))
        v1 = cheb_gauss_coeffs2vals(c[:, 1])
        v2 = cheb_gauss_coeffs2vals(c[:, 2])
        @test isapprox(v1, reverse(v1), atol=tol)
        @test isapprox(v2, -reverse(v2), atol=tol)

        v1_via_matrix = coeffs2vals_via_matrix(c[:, 1])
        v2_via_matrix = coeffs2vals_via_matrix(c[:, 2])
        @test isapprox(v1, reverse(v1), atol=tol)
        @test isapprox(v2, -reverse(v2), atol=tol)
    end

    @testset "Operator style" begin
        n = 100
        coeffs = rand(n)
        ctx = cheb_gauss_coeffs2vals_context(Float64, n)

        # Test operator call
        vals1 = cheb_gauss_coeffs2vals!(ctx, coeffs)
        vals2 = cheb_gauss_coeffs2vals(coeffs)
        @test isapprox(vals1, vals2, atol=tol)

        # Test multiple calls
        for _ in 1:10
            coeffs = rand(n)
            vals1 = cheb_gauss_coeffs2vals!(ctx, coeffs)
            vals2 = cheb_gauss_coeffs2vals(coeffs)
            @test isapprox(vals1, vals2, atol=tol)
        end
    end
end
