using Test
using GRSuite.ChebyshevSuite

@testset "Chebyshev Grid" begin
    @testset "First kind" begin
        @test cheb_grid(Float64, 5, FirstKind) ≈ [
            -0.9510565162951535,
            -0.5877852522924731,
            0.0,
            0.5877852522924731,
            0.9510565162951535,
        ]
    end

    @testset "Second kind" begin
        @test cheb_grid(Float64, 5, SecondKind) ≈
            [-1.0, -0.7071067811865475, 0.0, 0.7071067811865475, 1.0]
    end
end
