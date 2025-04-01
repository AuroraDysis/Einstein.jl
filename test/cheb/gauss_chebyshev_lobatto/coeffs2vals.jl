@testitem "gauss_chebyshev_lobatto_coeffs2vals" begin
    using LinearAlgebra

    # Set tolerance
    tol = typetol(Float64)

    # Test single coefficient conversion
    function coeffs2vals_via_matrix(c::AbstractVector{T}) where {T<:Number}
        n = length(c)
        A = gauss_chebyshev_lobatto_coeffs2vals_matrix(T, n)
        return A * c
    end

    @testset "Even case" begin
        c = collect(Float64, 1:5)
        vTrue = [3; -4 + sqrt(2); 3; -4 - sqrt(2); 15]
        v1 = gauss_chebyshev_lobatto_coeffs2vals(c)
        v2 = coeffs2vals_via_matrix(c)
        @test isapprox(v1, vTrue, atol=tol)
        @test isapprox(v2, vTrue, atol=tol)
    end

    @testset "Odd case" begin
        c = collect(Float64, 1:6)
        vTrue = [-3; 7 / 2; -(11 / 2) + sqrt(5); 7 / 2; -(11 / 2) - sqrt(5); 21]
        v1 = gauss_chebyshev_lobatto_coeffs2vals(c)
        v2 = coeffs2vals_via_matrix(c)
        @test isapprox(v1, vTrue, atol=tol)
        @test isapprox(v2, vTrue, atol=tol)
    end

    @testset "Symmetry preservation" begin
        c = kron(ones(10), Matrix{Float64}(I, 2, 2))
        c1 = @view c[:, 1]
        c2 = @view c[:, 2]
        v1 = gauss_chebyshev_lobatto_coeffs2vals(c1)
        v2 = gauss_chebyshev_lobatto_coeffs2vals(c2)
        @test isapprox(v1, reverse(v1), atol=tol)
        @test isapprox(v2, -reverse(v2), atol=tol)
    end

    @testset "Operator style" begin
        n = 100
        coeffs = rand(n)
        op = gauss_chebyshev_lobatto_coeffs2vals(Float64, n)

        # Test operator call
        vals1 = op(coeffs)
        vals2 = gauss_chebyshev_lobatto_coeffs2vals(coeffs)
        @test isapprox(vals1, vals2, atol=tol)

        # Test multiple calls
        for _ in 1:10
            coeffs = rand(n)
            vals1 = op(coeffs)
            vals2 = gauss_chebyshev_lobatto_coeffs2vals(coeffs)
            @test isapprox(vals1, vals2, atol=tol)
        end
    end
end
