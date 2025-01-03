@testset "cheb1_angles" begin
    @testset "n = 5" begin
        n = 5
        θ = cheb1_angles(n)

        @test length(θ) == n
        @test θ ≈ [
            2.82743338823081,
            2.19911485751286,
            1.57079632679490,
            0.942477796076938,
            0.314159265358979,
        ]
    end

    @testset "n = 6" begin
        n = 6
        θ = cheb1_angles(n)

        @test length(θ) == n
        @test θ ≈ [
            2.87979326579064,
            2.35619449019235,
            1.83259571459405,
            1.30899693899575,
            0.785398163397448,
            0.261799387799149,
        ]
    end
end
