@testitem "fdm_operator" begin
    using LinearAlgebra, BandedMatrices

    @testset "Operator Application" begin
        # Test on quadratic function f(x) = x²
        dx = 0.0001
        x = (-10 * dx):dx:(10 * dx)
        f = x .^ 2

        # First derivative should give 2x
        op1 = fdm_operator(1, 4, dx)
        df = op1 * f
        @test df ≈ 2 .* x

        # Second derivative should give 2
        op2 = fdm_operator(2, 4, dx)
        d2f = op2 * f
        @test isapprox(d2f, ones(length(x)) * 2, atol=1e-12)
    end

    @testset "Operator Properties" begin
        op = fdm_operator(2, 4, 0.1)

        # Test stencil properties
        @test length(op.weights) == 5  # 5-point stencil for 2nd der, 4th order
        @test sum(op.weights) < 10 * eps(Float64)  # sum of weights should be close to zero
        @test op.weights[1] == op.weights[end]  # symmetry
    end

    @testset "Matrix Representation" begin
        dx = 0.001
        x = (-10 * dx):dx:(10 * dx)
        n = length(x)
        op = fdm_operator(2, 4, dx)
        mat = fdm_operator_matrix(op, n)

        # Test matrix properties
        @test size(mat) == (n, n)
        @test isa(mat, BandedMatrix)

        # Test matrix application
        f = x .^ 2
        df = mat * f
        @test isapprox(df, 2 * ones(n), atol=1e-12)
    end
end

@testitem "Dissipation operator" begin
    using BandedMatrices

    @testset "Construction" begin
        op = fdm_dissipation_operator(2, 1.0, 0.1)
        @test op isa FiniteDifferenceOperator{Float64}
        @test op.factor[] ≈ 10.0
    end

    @testset "Weights correctness" begin
        op = fdm_dissipation_operator(4, 1.0, 1.0)
        expected_weights = fdm_dissipation_weights(4)
        @test all(op.weights .≈ expected_weights)
    end

    @testset "Matrix representation" begin
        n = 10
        op = fdm_dissipation_operator(fdm_dissipation_order(4), 1.0, 1.0)
        mat = fdm_operator_matrix(op, n)
        @test size(mat) == (n, n)
        @test isa(mat, BandedMatrix)
    end

    @testset "Array application" begin
        n = 10
        u = ones(n)
        op = fdm_dissipation_operator(fdm_dissipation_order(4), 1.0, 1.0)
        du = op * u
        @test all(abs.(du) .< 1e-12)
    end

    @testset "Smooth function" begin
        dx = 0.001
        u = sin.(0:dx:1.0)
        op = fdm_dissipation_operator(fdm_dissipation_order(4), 1.0, dx)
        du = op * u
        @test all(abs.(du) .< 1e-12)
    end
end
