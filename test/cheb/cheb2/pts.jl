using TestItems

@testitem "cheb2_pts" begin
    using GRSuite.ChebSuite, Test

    @testset "n = 5" begin
        n = 5
        x = cheb2_pts(n)

        @test length(x) == n
        @test x ≈ [-1.0, -0.707106781186548, 0.0, 0.707106781186548, 1.0]
    end

    @testset "n = 6" begin
        n = 6
        x = cheb2_pts(n)

        @test length(x) == n
        @test x ≈ [
            -1.0,
            -0.809016994374948,
            -0.309016994374947,
            0.309016994374947,
            0.809016994374948,
            1.0,
        ]
    end
end
