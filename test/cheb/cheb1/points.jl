using TestItems

@testitem "cheb1_points" begin
    using ApproxFun, Einstein.ChebyshevSuite, Test

    for TF in [Float64, BigFloat]
        tol = 10 * eps(TF)
        for intervals in [[-1, 1], [0, 1]]
            x_min = convert(TF, intervals[1])
            x_max = convert(TF, intervals[2])
            @testset "$TF, $intervals" begin
                for n in 0:10
                    @test isapprox(
                        cheb1_points(TF, n, x_min, x_max),
                        reverse(points(Chebyshev(x_min .. x_max), n)),
                        atol=tol,
                    )
                end
            end
        end
    end
end
