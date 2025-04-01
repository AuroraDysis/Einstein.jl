using TestItems

@testitem "gauss_chebyshev_lobatto_vals2coeffs" begin
    # Set tolerance
    tol = typetol(Float64)

    function vals2coeffs_via_matrix(v::AbstractVector{T}) where {T<:Number}
        n = length(v)
        A = gauss_chebyshev_lobatto_vals2coeffs_matrix(T, n)
        return A * v
    end

    @testset "Simple data" begin
        v = collect(Float64, 1:5)
        # Exact coefficients
        cTrue = [3.0, 1 + 1 / sqrt(2), 0.0, 1 - 1 / sqrt(2), 0.0]
        c1 = gauss_chebyshev_lobatto_vals2coeffs(v)
        c2 = vals2coeffs_via_matrix(v)
        @test isapprox(c1, cTrue, atol=tol)
        @test isapprox(c2, cTrue, atol=tol)
    end

    @testset "Array input" begin
        v = collect(Float64, 1:5)
        cTrue = [3.0, 1 + 1 / sqrt(2), 0.0, 1 - 1 / sqrt(2), 0.0]

        # Test forward and reversed arrays
        c1 = gauss_chebyshev_lobatto_vals2coeffs(v)
        c2 = gauss_chebyshev_lobatto_vals2coeffs(reverse(v))
        c3 = vals2coeffs_via_matrix(v)
        c4 = vals2coeffs_via_matrix(reverse(v))

        tmp = ones(length(cTrue))
        tmp[(end - 1):-2:1] .= -1

        @test isapprox(c1, cTrue, atol=tol)
        @test isapprox(c2, tmp .* cTrue, atol=tol)
        @test isapprox(c3, cTrue, atol=tol)
        @test isapprox(c4, tmp .* cTrue, atol=tol)
    end

    @testset "Symmetry preservation" begin
        # Create test data with even/odd symmetry
        n = 10
        v_even = repeat([1.0], n)
        v_odd = repeat([-1.0, 1.0], n ÷ 2)

        c_even = gauss_chebyshev_lobatto_vals2coeffs(v_even)
        c_odd = gauss_chebyshev_lobatto_vals2coeffs(v_odd)

        # Even coefficients should have zero odd terms
        @test all(abs.(c_even[2:2:end]) .< tol)
        # Odd coefficients should have zero even terms
        @test all(abs.(c_odd[1:2:end]) .< tol)
    end

    @testset "Operator style" begin
        n = 100
        vals = rand(n)
        op = GaussChebyshevLobattoGrid.vals2coeffs(Float64, n)

        # Test operator call
        coeffs1 = op(vals)
        coeffs2 = gauss_chebyshev_lobatto_vals2coeffs(vals)
        @test isapprox(coeffs1, coeffs2, atol=tol)

        # Test multiple calls
        for _ in 1:10
            vals = rand(n)
            coeffs1 = op(vals)
            coeffs2 = gauss_chebyshev_lobatto_vals2coeffs(vals)
            @test isapprox(coeffs1, coeffs2, atol=tol)
        end
    end
end