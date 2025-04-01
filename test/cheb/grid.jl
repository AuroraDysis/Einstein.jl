using TestItems

@testitem "ChebyshevGrid" begin
    using ApproxFun, Einstein.ChebyshevSuite, Test

    for TF in [Float64, BigFloat]
        lower_bound = zero(TF)
        upper_bound = one(TF)
        @testset "$TF, FirstKind" begin
            for n in 0:10
                grid = ChebyshevGrid(n, lower_bound, upper_bound; kind=1)
                @test isapprox(
                    grid.data, chebtech1_points(TF, n, lower_bound, upper_bound), atol=10 * eps(TF)
                )
            end
        end

        @testset "$TF, SecondKind" begin
            for n in 0:10
                grid = ChebyshevGrid(n, lower_bound, upper_bound; kind=2)
                @test isapprox(
                    grid.data, chebtech2_points(TF, n, lower_bound, upper_bound), atol=10 * eps(TF)
                )
            end
        end
    end
end
