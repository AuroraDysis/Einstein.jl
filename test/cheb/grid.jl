using TestItems

@testitem "cheb_grid" begin
    using ApproxFun, Einstein.ChebSuite, Test

    for TR in [Float64, BigFloat]
        x_min = zero(TR)
        x_max = one(TR)
        @testset "$TR" begin
            for n in 0:10
                grid = cheb_grid(x_min, x_max, n; node=ChebyshevNode.FirstKind)
                @test isapprox(grid.x, cheb1_pts(TR, n, x_min, x_max), atol=10 * eps(TR))
            end
            for n in 0:10
                grid = cheb_grid(x_min, x_max, n; node=ChebyshevNode.SecondKind)
                @test isapprox(grid.x, cheb2_pts(TR, n, x_min, x_max), atol=10 * eps(TR))
            end
        end
    end
end
