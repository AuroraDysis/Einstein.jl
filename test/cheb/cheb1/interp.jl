using TestItems

@testitem "cheb1_interp" begin
    using PDESuite.ChebSuite, Test

    tol = 100 * eps()

    n = 5
    x = cheb1_pts(n)
    v = sin.(x)
    op = Cheb1InterpOp(n)

    for i in 1:n
        @test op(v, x[i]) ≈ v[i]
    end
end

@testitem "Interpolation accuracy" begin
    using PDESuite.ChebSuite, Test

    tol = 100 * eps()

    n = 32
    op = Cheb1InterpOp(n)
    x = cheb1_pts(n)
    v = cos.(x)

    # Test at intermediate points
    test_points = [-0.9, -0.5, 0.0, 0.5, 0.9]
    for xi in test_points
        y_interp = op(v, xi)
        y_exact = cos(xi)
        @test abs(y_interp - y_exact) < tol
    end
end
