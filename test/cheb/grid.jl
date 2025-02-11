using TestItems

@testitem "ChebyshevGrid" begin
    using ApproxFun, Einstein.ChebyshevSuite, Test

    for TR in [Float64, BigFloat]
        x_min = zero(TR)
        x_max = one(TR)
        @testset "$TR, FirstKind" begin
            for n in 0:10
                grid = ChebyshevGrid(n, x_min, x_max, ChebyshevNode.FirstKind)
                @test isapprox(grid.data, cheb1_points(TR, n, x_min, x_max), atol=10 * eps(TR))
            end
        end

        @testset "$TR, SecondKind" begin
            for n in 0:10
                grid = ChebyshevGrid(n, x_min, x_max, ChebyshevNode.SecondKind)
                @test isapprox(grid.data, cheb2_points(TR, n, x_min, x_max), atol=10 * eps(TR))
            end
        end
    end
end
