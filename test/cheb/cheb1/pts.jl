@testset "cheb1_pts" begin
    for TR in [Float64, BigFloat]
        tol = 10 * eps(TR)
        for intervals in [[-1, 1], [0, 1]]
            x_min = convert(TR, intervals[1])
            x_max = convert(TR, intervals[2])
            @testset "$TR, $intervals" begin
                for n in 0:10
                    @test isapprox(
                        cheb1_pts(TR, n, x_min, x_max),
                        reverse(points(Chebyshev(x_min .. x_max), n)),
                        atol=tol,
                    )
                end
            end
        end
    end
end
