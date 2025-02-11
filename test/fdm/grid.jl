@testitem "fdm_grid" begin
    using Einstein.FiniteDifferenceSuite, Test

    @testset "Basic functionality" begin
        # Integer-spaced grid
        grid = UniformGrid(0.0, 1.0, 0.1)
        @test length(grid) == 11
        @test grid[1] ≈ 0.0
        @test grid[end] ≈ 1.0
        @test all(diff(grid) .≈ 0.1)
    end

    @testset "Edge cases" begin
        # Small spacing
        grid = UniformGrid(0.0, 1.0, 1e-5)
        @test length(grid) == 100001
        @test grid[1] == 0.0
        @test isapprox(grid[end], 1.0, atol=10 * eps(Float64))

        # Large interval
        grid = UniformGrid(0.0, 1000.0, 0.1)
        @test length(grid) == 10001
        @test grid[1] == 0.0
        @test isapprox(grid[end], 1000.0, atol=10 * eps(Float64))

        # Float32 type
        grid = UniformGrid(0.0f0, 1.0f0, 0.1f0)
        @test eltype(grid) == Float32
        @test isapprox(grid[end], 1.0f0, atol=10 * eps(Float32))
    end

    @testset "Error cases" begin
        # Invalid interval
        @test_throws ArgumentError UniformGrid(1.0, 0.0, 0.1)

        # Invalid spacing
        @test_throws ArgumentError UniformGrid(0.0, 1.0, 0.0)
        @test_throws ArgumentError UniformGrid(0.0, 1.0, -0.1)

        # Test endpoint mismatch case
        # This should throw an error when dx doesn't divide interval evenly
        @test_throws ArgumentError UniformGrid(0.0, 1.0, 0.3)
    end
end