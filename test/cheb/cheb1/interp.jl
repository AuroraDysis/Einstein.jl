@testitem "cheb1_interp" begin
    using Einstein.ChebSuite, Test

    @testset "real" begin
        tol = 100 * eps()
        n = 40
        x = cheb1_pts(n)
        v = sin.(x)
        op = Cheb1InterpOp(n)
    
        for i in 1:n
            @test op(v, x[i]) ≈ v[i]
        end

        for i in 1:n-1
            x0 = (x[i] + x[i+1]) / 2
            @test isapprox(op(v, x0), sin(x0), atol=tol)
        end
    end

    @testset "complex" begin
        tol = 100 * eps()
        n = 30
        x = cheb1_pts(n)
        v = sin.(x) + im * cos.(x)
        op = Cheb1InterpOp(n)
    
        for i in 1:n
            @test op(v, x[i]) ≈ v[i]
        end

        for i in 1:n-1
            x0 = (x[i] + x[i+1]) / 2
            @test isapprox(op(v, x0), sin(x0) + im * cos(x0), atol=tol)
        end
    end
end
