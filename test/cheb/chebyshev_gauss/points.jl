@testitem "cheb_gauss_points" begin
    for TF in [Float64, BigFloat]
        tol = 10 * eps(TF)
        for intervals in [[-1, 1], [0, 1]]
            lower_bound = convert(TF, intervals[1])
            upper_bound = convert(TF, intervals[2])
            @testset "$TF, $intervals" begin
                for n in 0:10
                    x = cheb_gauss_points(TF, n, lower_bound, upper_bound)
                    
                    # Test number of points
                    @test length(x) == max(0, n)
                    
                    # Test bounds
                    if n > 0
                        @test minimum(x) >= lower_bound
                        @test maximum(x) <= upper_bound
                    end
                    
                    # Test symmetry for [-1,1] interval
                    if intervals == [-1, 1] && n > 1
                        @test isapprox(x, -reverse(x), atol=tol)
                    end
                    
                    # Test spacing between points
                    if n > 1
                        dx = diff(x)
                        @test all(dx .> 0)  # Points should be strictly increasing
                    end
                end
            end
        end
    end
end
