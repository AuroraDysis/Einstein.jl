@testitem "fdm_operator" begin
    using LinearAlgebra, Einstein.FiniteDifferenceSuite, Test

    @testset "Operator Application" begin
        # Test on quadratic function f(x) = x²
        dx = 0.001
        x = (-10 * dx):dx:(10 * dx)
        f = x .^ 2
        center_idx = length(x) ÷ 2

        # First derivative should give 2x
        op1 = fdm_operator(Float64, 1, 4, dx)
        df = op1 * f
        @test df ≈ 2 .* x

        # Second derivative should give 2
        op2 = fdm_operator(Float64, 2, 4, dx)
        d2f = op2 * f
        @test isapprox(d2f[center_idx], 2, rtol=1e-14)
    end

    @testset "Operator Properties" begin
        op = fdm_operator(Float64, 2, 4, 0.1)

        # Test stencil properties
        @test length(op.weights) == 5  # 5-point stencil for 2nd der, 4th order
        @test sum(op.weights) < 10 * eps(Float64)  # sum of weights should be close to zero
        @test op.weights[1] == op.weights[end]  # symmetry
    end
end
