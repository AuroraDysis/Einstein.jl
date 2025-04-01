@testitem "GaussChebyshevGrid - points" begin
    using ApproxFun

    for TF in [Float64, BigFloat]
        tol = 10 * eps(TF)
        for intervals in [[-1, 1], [0, 1]]
            lower_bound = convert(TF, intervals[1])
            upper_bound = convert(TF, intervals[2])
            @testset "$TF, $intervals" begin
                for n in 0:10
                    @test isapprox(
                        GaussChebyshevGrid.points(TF, n, lower_bound, upper_bound),
                        reverse(points(Chebyshev(lower_bound .. upper_bound), n)),
                        atol=tol,
                    )
                end
            end
        end
    end
end
