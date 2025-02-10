using TestItems

@testitem "cheb_grid" begin
    using ApproxFun, Einstein.ChebSuite, Test

    for TR in [Float64, BigFloat]
        x_min = zero(TR)
        x_max = one(TR)
        for TNode in [ChebyshevNode.FirstKind, ChebyshevNode.SecondKind]
            @testset "$TR, $TNode" begin
                for n in 0:10
                    grid = cheb_grid(TNode, n, x_min, x_max)
                    @test isapprox(grid.data, cheb_pts(TR, TNode, n, x_min, x_max), atol=10 * eps(TR))
                end
            end
        end
    end
end
