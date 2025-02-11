using TestItems

@testitem "ChebyshevSynthesis and ChebyshevAnalysis" begin
    using LinearAlgebra, Einstein.ChebyshevSuite, Test

    for kind in [1, 2]
        @testset "kind = $kind, real" begin
            tol = 100 * eps()
            n = 40

            grid = ChebyshevGrid(n, 0.0, 1.0; kind=kind)
            syn = ChebyshevSynthesis(grid)
            ana = ChebyshevAnalysis(grid)
            v = sin.(grid)

            c = ana(v)
            v1 = syn(c)
            @test norm(v - v1, Inf) < tol
        end

        @testset "kind = $kind, complex" begin
            tol = 100 * eps()
            n = 40
            grid = ChebyshevGrid(n, 0.0, 1.0; kind=kind)
            syn = ChebyshevSynthesis(grid)
            ana = ChebyshevAnalysis(grid)
            v = sin.(grid) + im * cos.(grid)

            c = ana(v)
            v1 = syn(c)
            @test norm(v - v1, Inf) < tol
        end
    end
end
