@testitem "FDM Convolution Operations" begin
    using StaticArrays

    @testset "Interior Convolution" begin
        # Test 1D interior convolution with constant factor
        in_1d = [1.0, 2.0, 3.0, 4.0, 5.0]
        out_1d = similar(in_1d)
        weights = SVector{3}([-0.5, 0.0, 0.5])
        factor = 1.0

        fdm_convolve_interior!(out_1d, in_1d, weights, factor)
        @test out_1d[2:4] ≈ [1.0, 1.0, 1.0]

        # Test 2D interior convolution
        in_2d = [1.0 2.0; 2.0 3.0; 3.0 4.0; 4.0 5.0]
        out_2d = similar(in_2d)
        weights = SVector{3}([-0.5, 0.0, 0.5])
        
        fdm_convolve_interior!(out_2d, in_2d, weights, factor)
        @test out_2d[2:3, :] ≈ [1.0 1.0; 1.0 1.0]
    end

    @testset "Boundary Convolution" begin
        # Test 1D boundary convolution
        in_1d = [1.0, 2.0, 3.0, 4.0, 5.0]
        out_1d = similar(in_1d)
        left_weights = @SMatrix [1.0 -2.0; -0.5 0.5]
        right_weights = @SMatrix [-0.5 0.5; 1.0 -2.0]
        factor = 1.0

        fdm_convolve_boundary!(out_1d, in_1d, left_weights, right_weights, factor)
        @test out_1d[1:2] ≈ [-3.0, 0.5]
        @test out_1d[4:5] ≈ [0.5, -6.0]

        # Test with vector factors
        factors = [0.5, 1.0, 1.0, 1.0, 0.5]
        fdm_convolve_boundary!(out_1d, in_1d, left_weights, right_weights, factors)
        @test out_1d[1:2] ≈ [-1.5, 0.5]
        @test out_1d[4:5] ≈ [0.5, -3.0]
    end
end
