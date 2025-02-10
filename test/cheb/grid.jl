using TestItems

@testitem "cheb_grid" begin
    using ApproxFun, Einstein.ChebSuite, Test

    for TR in [Float64, BigFloat]
        x_min = zero(TR)
        x_max = one(TR)
        @testset "$TR, ChebNode.FirstKind" begin
            for n in 0:10
                grid = cheb_grid(ChebNode.FirstKind, n, x_min, x_max)
                @test isapprox(grid.data, cheb1_pts(TR, n, x_min, x_max), atol=10 * eps(TR))
            end
        end

        @testset "$TR, ChebNode.SecondKind" begin
            for n in 0:10
                grid = cheb_grid(ChebNode.SecondKind, n, x_min, x_max)
                @test isapprox(grid.data, cheb2_pts(TR, n, x_min, x_max), atol=10 * eps(TR))
            end
        end
    end
end
