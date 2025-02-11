using TestItems

@testitem "ChebyshevGrid" begin
    using ApproxFun, Einstein.ChebyshevSuite, Test

    for TF in [Float64, BigFloat]
        x_min = zero(TF)
        x_max = one(TF)
        @testset "$TF, FirstKind" begin
            for n in 0:10
                grid = ChebyshevGrid(n, x_min, x_max, ChebyshevNode.FirstKind)
                @test isapprox(grid.data, cheb1_points(TF, n, x_min, x_max), atol=10 * eps(TF))
            end
        end

        @testset "$TF, SecondKind" begin
            for n in 0:10
                grid = ChebyshevGrid(n, x_min, x_max, ChebyshevNode.SecondKind)
                @test isapprox(grid.data, cheb2_points(TF, n, x_min, x_max), atol=10 * eps(TF))
            end
        end
    end
end
