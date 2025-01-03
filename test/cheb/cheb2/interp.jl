using TestItems

@testitem "cheb2_interp" begin
    using GRSuite.ChebSuite, Test

    tol = 100 * eps()
    n = 5
    x = cheb2_pts(n)
    v = sin.(x)
    op = Cheb2InterpOp(n)

    for i in 1:n
        @test op(v, x[i]) ≈ v[i]
    end
end
