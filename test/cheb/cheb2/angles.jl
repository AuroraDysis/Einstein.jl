using TestItems

@testitem "chebtech2_angles" begin
    using Einstein.ChebyshevSuite, Test

    @testset "n = 5" begin
        n = 5
        θ = chebtech2_angles(n)

        @test length(θ) == n
        @test θ ≈
            [3.14159265358979, 2.35619449019235, 1.57079632679490, 0.785398163397448, 0.0]
    end

    @testset "n = 6" begin
        n = 6
        θ = chebtech2_angles(n)

        @test length(θ) == n
        @test θ ≈ [
            3.14159265358979,
            2.51327412287183,
            1.88495559215388,
            1.25663706143592,
            0.628318530717959,
            0.0,
        ]
    end
end
