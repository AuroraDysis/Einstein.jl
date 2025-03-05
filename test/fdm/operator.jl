@testitem "fdm_derivative_operator" begin
    using LinearAlgebra, BandedMatrices

    @testset "Operator Application" begin
        # Test on quadratic function f(x) = x²
        dx = 0.0001
        x = (-10 * dx):dx:(10 * dx)
        f = x .^ 2

        # First derivative should give 2x
        op1 = fdm_derivative_operator(1, 4, dx)
        df = op1 * f
        @test df ≈ 2 .* x

        # Second derivative should give 2
        op2 = fdm_derivative_operator(2, 4, dx)
        d2f = op2 * f
        @test isapprox(d2f, ones(length(x)) * 2, atol=1e-12)
    end

    @testset "Operator Properties" begin
        op = fdm_derivative_operator(2, 4, 0.1)

        # Test stencil properties
        @test length(op.weights) == 5  # 5-point stencil for 2nd der, 4th order
        @test sum(op.weights) < 10 * eps(Float64)  # sum of weights should be close to zero
        @test op.weights[1] == op.weights[end]  # symmetry
    end

    @testset "Matrix Representation" begin
        dx = 0.001
        x = (-10 * dx):dx:(10 * dx)
        n = length(x)
        op = fdm_derivative_operator(2, 4, dx)
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

@testitem "fdm_hermite_derivative_operator" begin
    @testset "Operator Construction" begin
        # Test basic construction
        dx = 0.1
        op = fdm_hermite_derivative_operator(2, 4, dx)
        
        # Test operator properties
        @test op isa HermiteFiniteDifferenceOperator{Float64}
        @test op.derivative_order == 2
        @test op.accuracy_order == 4
        
        # Check that weights are consistent with fdm_hermite_weights
        D_weights, E_weights = fdm_hermite_weights(2, 4)
        @test all(op.D_weights .≈ D_weights)
        @test all(op.E_weights .≈ E_weights)
    end
    
    @testset "Operator Application" begin
        # Test on a cubic function f(x) = x³
        # For f(x) = x³, f'(x) = 3x², f''(x) = 6x
        dx = 0.001
        x = (-5 * dx):dx:(5 * dx)
        f = x .^ 3
        df = 3 .* x .^ 2
        
        # Create the operator for second derivative
        op = fdm_hermite_derivative_operator(2, 4, dx)
        
        # Apply the operator with both function and derivative
        # To use the Hermite operator, we need to provide both f and df
        # Since the operator doesn't have a direct multiplication method with two inputs,
        # we'll create temporary arrays and use the apply function
        d2f = similar(f)
        fdm_apply_operator!(d2f, op, f, df)
        
        # Second derivative should be 6x
        @test d2f ≈ 6 .* x atol=1e-10
    end
end

@testitem "Dissipation operator" begin
    using BandedMatrices

    @testset "Construction" begin
        op = fdm_dissipation_operator(2, 1.0, 0.1)
        @test op isa FiniteDifferenceDissipationOperator{Float64}
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
