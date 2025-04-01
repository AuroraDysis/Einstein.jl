using TestItems

@testitem "GaussChebyshevLobatto - vals2coeffs" begin
    # Set tolerance
    tol = typetol(Float64)

    function vals2coeffs_via_matrix(v::AbstractVector{T}) where {T<:Number}
        n = length(v)
        A = GaussChebyshevLobatto.vals2coeffs_matrix(T, n)
        return A * v
    end

    @testset "Single value" begin
        v = sqrt(2)
        c1 = GaussChebyshevLobatto.vals2coeffs([v])
        c2 = vals2coeffs_via_matrix([v])
        @test v ≈ c1[1]
        @test v ≈ c2[1]
    end

    @testset "Simple data" begin
        v = collect(1.0:5.0)
        # Exact coefficients
        cTrue = [3.0, 1 + 1 / sqrt(2), 0.0, 1 - 1 / sqrt(2), 0.0]
        c1 = GaussChebyshevLobatto.vals2coeffs(v)
        c2 = vals2coeffs_via_matrix(v)
        @test isapprox(c1, cTrue, atol=tol)
        @test isapprox(c2, cTrue, atol=tol)
    end

    @testset "Array input" begin
        v = collect(1.0:5.0)
        cTrue = [3.0, 1 + 1 / sqrt(2), 0.0, 1 - 1 / sqrt(2), 0.0]

        # Test forward and reversed arrays
        c1 = GaussChebyshevLobatto.vals2coeffs(v)
        c2 = GaussChebyshevLobatto.vals2coeffs(reverse(v))
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

        c_even = GaussChebyshevLobatto.vals2coeffs(v_even)
        c_odd = GaussChebyshevLobatto.vals2coeffs(v_odd)

        # Even coefficients should have zero odd terms
        @test all(abs.(c_even[2:2:end]) .< tol)
        # Odd coefficients should have zero even terms
        @test all(abs.(c_odd[1:2:end]) .< tol)
    end

    @testset "Operator style" begin
        n = 100
        vals = rand(n)
        op = GaussChebyshevLobatto.vals2coeffs(Float64, n)

        # Test operator call
        coeffs1 = op(vals)
        coeffs2 = GaussChebyshevLobatto.vals2coeffs(vals)
        @test isapprox(coeffs1, coeffs2, atol=tol)

        # Test multiple calls
        for _ in 1:10
            vals = rand(n)
            coeffs1 = op(vals)
            coeffs2 = GaussChebyshevLobatto.vals2coeffs(vals)
            @test isapprox(coeffs1, coeffs2, atol=tol)
        end
    end
end