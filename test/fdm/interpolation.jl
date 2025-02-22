@testitem "LocalBarycentricInterpolation" begin
    @testset "Constructor validation" begin
        points = 0.0:0.1:1.0
        values = sin.(points)
        
        # Valid construction
        @test_nowarn LocalBarycentricInterpolation(points, values, degree=4)
        
        # Invalid degree (too few points)
        too_few_points = 0.0:0.1:0.2  # only 3 points
        @test_throws ArgumentError LocalBarycentricInterpolation(too_few_points, sin.(too_few_points), degree=4)
        
        # Mismatched lengths
        @test_throws ArgumentError LocalBarycentricInterpolation(points, values[1:end-1], degree=4)
    end

    @testset "Interpolation accuracy" begin
        # Test exact interpolation at points
        points = 0.0:0.01:1.0
        f(x) = sin(π*x)
        values = f.(points)
        itp = LocalBarycentricInterpolation(points, values, degree=4)
        
        for x in points
            @test isapprox(itp(x), f(x), rtol=1e-14)
        end
        
        # Test interpolation between points
        x_test = 0.05:0.1:0.95
        for x in x_test
            @test isapprox(itp(x), f(x), rtol=1e-14)
        end
    end

    @testset "Different polynomial degrees" begin
        points = 0.0:0.01:1.0
        f(x) = x^2  # quadratic function
        values = f.(points)
        
        for degree in [2, 3, 4, 5]
            itp = LocalBarycentricInterpolation(points, values, degree=degree)
            x_test = 0.05
            @test isapprox(itp(x_test), f(x_test), rtol=1e-14)
        end
    end

    @testset "Boundary behavior" begin
        points = 0.0:0.1:1.0
        values = sin.(points)
        itp = LocalBarycentricInterpolation(points, values, degree=4)
        
        # Test at boundaries
        @test isapprox(itp(first(points)), first(values), rtol=1e-14)
        @test isapprox(itp(last(points)), last(values), rtol=1e-14)
        
        # Test near boundaries
        @test_nowarn itp(0.01)
        @test_nowarn itp(0.99)
    end

    @testset "Simple functions" begin
        points = 0.0:0.01:1.0
        
        # Linear function
        f_linear(x) = 2x + 1
        values = f_linear.(points)
        itp = LocalBarycentricInterpolation(points, values, degree=2)
        @test isapprox(itp(0.55), f_linear(0.55), rtol=1e-14)
        
        # Quadratic function
        f_quad(x) = x^2 - x + 1
        values = f_quad.(points)
        itp = LocalBarycentricInterpolation(points, values, degree=3)
        @test isapprox(itp(0.55), f_quad(0.55), rtol=1e-14)
    end
end