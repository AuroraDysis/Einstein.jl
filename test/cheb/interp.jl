@testitem "cheb1_interp" begin
    using Einstein.ChebyshevSuite, Test

    for kind in [ChebyshevFirstKindNode, ChebyshevSecondKindNode]
        @testset "kind = $kind, real" begin
            tol = 100 * eps()
            n = 40
    
            grid = cheb_grid(kind, n, -1.0, 1.0)
            op = cheb_interp(grid)
            v = sin.(grid)
        
            for i in 1:n
                @test op(v, grid[i]) ≈ v[i]
            end
    
            for i in 1:n-1
                x0 = (grid[i] + grid[i+1]) / 2
                @test isapprox(op(v, x0), sin(x0), atol=tol)
            end
        end
    
        @testset "kind = $kind, complex" begin
            tol = 100 * eps()
            n = 30
            grid = cheb_grid(kind, n, -1.0, 1.0)
            op = cheb_interp(grid)
            v = sin.(grid) + im * cos.(grid)
        
            for i in 1:n
                @test op(v, grid[i]) ≈ v[i]
            end
    
            for i in 1:n-1
                x0 = (grid[i] + grid[i+1]) / 2
                @test isapprox(op(v, x0), sin(x0) + im * cos(x0), atol=tol)
            end
        end
    end
end
