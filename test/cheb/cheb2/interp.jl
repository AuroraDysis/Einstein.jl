@testset "cheb2_interp" begin
    tol = 100 * eps()

    @testset "Interpolation at nodes" begin
        n = 5
        x = cheb2_pts(n)
        v = sin.(x)
        op = Cheb2InterpOp(n)

        for i in 1:n
            @test op(v, x[i]) ≈ v[i]
        end
    end

    @testset "Interpolation accuracy" begin
        n = 32
        op = Cheb2InterpOp(n)
        x = cheb2_pts(n)
        v = cos.(x)

        # Test at intermediate points
        test_points = [-0.9, -0.5, 0.0, 0.5, 0.9]
        for xi in test_points
            y_interp = op(v, xi)
            y_exact = cos(xi)
            @test abs(y_interp - y_exact) < tol
        end
    end
end
