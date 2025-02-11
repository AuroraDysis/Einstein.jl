@testitem "ChebyshevGrid BarycentricInterpolation" begin
    using Einstein.ChebyshevSuite, Test

    for kind in [1, 2]
        @testset "kind = $kind, real" begin
            tol = 100 * eps()
            n = 40

            grid = ChebyshevGrid(n, -1.0, 1.0; kind=kind)
            itp = ChebyshevInterpolation(grid)
            v = sin.(grid)

            for i in 1:n
                @test itp(v, grid[i]) ≈ v[i]
            end

            for i in 1:(n - 1)
                x0 = (grid[i] + grid[i + 1]) / 2
                @test isapprox(itp(v, x0), sin(x0), atol=tol)
            end
        end

        @testset "kind = $kind, complex" begin
            tol = 100 * eps()
            n = 30
            grid = ChebyshevGrid(n, -1.0, 1.0; kind=kind)
            itp = ChebyshevInterpolation(grid)
            v = sin.(grid) + im * cos.(grid)

            for i in 1:n
                @test itp(v, grid[i]) ≈ v[i]
            end

            for i in 1:(n - 1)
                x0 = (grid[i] + grid[i + 1]) / 2
                @test isapprox(itp(v, x0), sin(x0) + im * cos(x0), atol=tol)
            end
        end
    end
end
