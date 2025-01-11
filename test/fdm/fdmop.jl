@testitem "fdm_centralop" begin
    using LinearAlgebra, Einstein.FDMSuite, Test

    @testset "Constructor Tests" begin
        dx = 0.1
        # Basic construction
        op = fdm_centralop(2, 4, dx)  # 2nd derivative, 4th order accuracy
        @test op.der_order == 2
        @test op.acc_order == 4
    end

    @testset "Operator Application" begin
        # Test on quadratic function f(x) = x²
        dx = 0.001
        x = -10*dx:dx:10*dx
        f = x .^ 2
        center_idx = length(x) ÷ 2

        # First derivative should give 2x
        op1 = fdm_centralop(1, 4, dx)
        df = op1(f[center_idx-op1.num_side:center_idx+op1.num_side])
        @test df ≈ 2 .* x[center_idx]

        # Second derivative should give 2
        op2 = fdm_centralop(2, 4, dx)
        d2f = op2(f[center_idx-op2.num_side:center_idx+op2.num_side])
        @test d2f ≈ 2.0
    end

    @testset "Operator Properties" begin
        op = fdm_centralop(2, 4, 0.1)

        # Test stencil properties
        @test length(op.wts) == 5  # 5-point stencil for 2nd der, 4th order
        @test sum(op.wts) < 10 * eps(Float64)  # sum of weights should be close to zero
        @test op.wts[1] == op.wts[end]  # symmetry
    end
end

@testitem "fdm_hermiteop" begin
    using LinearAlgebra, Einstein.FDMSuite, Test

    @testset "Constructor Tests" begin
        dx = 0.1
        # Basic construction
        op = fdm_hermiteop(2, 4, dx)  # 2nd derivative, 4th order accuracy
        @test op.der_order == 2
        @test op.acc_order == 4
    end

    @testset "Operator Application" begin
        # Test on quadratic function f(x) = x²
        dx = 0.001
        x = collect(-10*dx:dx:10*dx)
        f = x .^ 2
        df = 2 .* x
        center_idx = length(x) ÷ 2

        # Second derivative should give 2
        op2 = fdm_hermiteop(2, 4, dx)
        d2f = op2(f[center_idx-op2.num_side:center_idx+op2.num_side], df[center_idx-op2.num_side:center_idx+op2.num_side])
        @test d2f ≈ 2.0
    end
end
