@testitem "ChebyshevSuite - BarycentricInterpolation" begin
    for grid in [GaussChebyshevGrid, GaussChebyshevLobattoGrid]
        @testset "grid = $grid, real" begin
            tol = typetol(Float64)
            n = 40

            points = grid.points(Float64, n, -1.0, 1.0)
            weights = grid.barycentric_weights(n)
            itp = BarycentricInterpolation(points, weights)
            v = sin.(points)

            for i in 1:n
                @test itp(v, points[i]) ≈ v[i]
            end

            for i in 1:(n - 1)
                x0 = (points[i] + points[i + 1]) / 2
                @test isapprox(itp(v, x0), sin(x0), atol=tol)
            end
        end

        @testset "complex" begin
            tol = 100 * eps()
            n = 30
            points = grid.points(Float64, n, -1.0, 1.0)
            weights = grid.barycentric_weights(n)
            itp = BarycentricInterpolation(points, weights)
            v = sin.(points) + im * cos.(points)

            for i in 1:n
                @test itp(v, points[i]) ≈ v[i]
            end

            for i in 1:(n - 1)
                x0 = (points[i] + points[i + 1]) / 2
                @test isapprox(itp(v, x0), sin(x0) + im * cos(x0), atol=tol)
            end
        end
    end
end
