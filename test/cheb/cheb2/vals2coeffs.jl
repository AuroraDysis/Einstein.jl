using TestItems

@testitem "cheb2_vals2coeffs" begin
    using GRSuite.ChebSuite, Test

    # Set tolerance
    tol = 100 * eps()

    @testset "Single value" begin
        v = sqrt(2)
        c = cheb2_vals2coeffs([v])
        @test v ≈ c[1]
    end

    @testset "Simple data" begin
        v = collect(1.0:5.0)
        # Exact coefficients
        cTrue = [3.0, 1 + 1 / sqrt(2), 0.0, 1 - 1 / sqrt(2), 0.0]
        c = cheb2_vals2coeffs(v)
        @test maximum(abs.(c .- cTrue)) < tol
        @test all(abs.(imag.(c)) .== 0)
    end

    @testset "Array input" begin
        v = collect(1.0:5.0)
        cTrue = [3.0, 1 + 1 / sqrt(2), 0.0, 1 - 1 / sqrt(2), 0.0]

        # Test forward and reversed arrays
        c1 = cheb2_vals2coeffs(v)
        c2 = cheb2_vals2coeffs(reverse(v))

        tmp = ones(length(cTrue))
        tmp[(end - 1):-2:1] .= -1

        @test maximum(abs.(c1 .- cTrue)) < tol
        @test maximum(abs.(c2 .- (tmp .* cTrue))) < tol
    end

    @testset "Symmetry preservation" begin
        # Create test data with even/odd symmetry
        n = 10
        v_even = repeat([1.0], n)
        v_odd = repeat([-1.0, 1.0], n ÷ 2)

        c_even = cheb2_vals2coeffs(v_even)
        c_odd = cheb2_vals2coeffs(v_odd)

        # Even coefficients should have zero odd terms
        @test all(abs.(c_even[2:2:end]) .< tol)
        # Odd coefficients should have zero even terms
        @test all(abs.(c_odd[1:2:end]) .< tol)
    end

    @testset "Operator style" begin
        n = 100
        vals = rand(n)
        op = Cheb2Vals2CoeffsOp(Float64, n)

        # Test operator call
        coeffs1 = op(vals)
        coeffs2 = cheb2_vals2coeffs(vals)
        @test maximum(abs.(coeffs1 .- coeffs2)) < tol

        # Test multiple calls
        for _ in 1:10
            vals = rand(n)
            coeffs1 = op(vals)
            coeffs2 = cheb2_vals2coeffs(vals)
            @test maximum(abs.(coeffs1 .- coeffs2)) < tol
        end
    end
end