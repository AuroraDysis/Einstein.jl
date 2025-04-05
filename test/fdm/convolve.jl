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
        in_2d = [1.0 2.0; 2.0 3.0; 3.0 4.0; 4.0 5.0; 5.0 6.0]
        out_2d = similar(in_2d)
        weights = SVector{3}([-0.5, 0.0, 0.5])
        
        fdm_convolve_interior!(out_2d, in_2d, weights, factor)
        @test out_2d[2:4, :] ≈ [1.0 1.0; 1.0 1.0; 1.0 1.0]

        # Test Add mode for 1D interior convolution
        out_1d_add = ones(size(in_1d))
        fdm_convolve_interior!(out_1d_add, in_1d, weights, factor, ConvolveAdd())
        @test out_1d_add[2:4] ≈ [2.0, 2.0, 2.0]  # original ones + convolution result

        # Test Add mode for 2D interior convolution
        out_2d_add = ones(size(in_2d))
        fdm_convolve_interior!(out_2d_add, in_2d, weights, factor, ConvolveAdd())
        @test out_2d_add[2:4, :] ≈ [2.0 2.0; 2.0 2.0; 2.0 2.0]
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

        # Test 2D boundary convolution
        in_2d = [1.0 2.0; 2.0 3.0; 3.0 4.0; 4.0 5.0; 5.0 6.0; 6.0 7.0]
        out_2d = similar(in_2d)
        left_weights = @SMatrix [1.0 -2.0; -0.5 0.5]
        right_weights = @SMatrix [-0.5 0.5; 1.0 -2.0]

        fdm_convolve_boundary!(out_2d, in_2d, left_weights, right_weights, factor)
        @test out_2d[1:2, :] ≈ [-3.0 -4.0; 0.5 0.5]
        @test out_2d[5:6, :] ≈ [0.5 0.5; -7.0 -8.0]

        # Test Add mode for 1D boundary convolution
        out_1d_add = ones(size(in_1d))
        fdm_convolve_boundary!(out_1d_add, in_1d, left_weights, right_weights, factor, ConvolveAdd())
        @test out_1d_add[1:2] ≈ [-2.0, 1.5]  # original ones + convolution result
        @test out_1d_add[4:5] ≈ [1.5, -5.0]

        # Test Add mode for 2D boundary convolution
        out_2d_add = ones(size(in_2d))
        fdm_convolve_boundary!(out_2d_add, in_2d, left_weights, right_weights, factor, ConvolveAdd())
        @test out_2d_add[1:2, :] ≈ [-2.0 -3.0; 1.5 1.5]
        @test out_2d_add[5:6, :] ≈ [1.5 1.5; -6.0 -7.0]
    end
end
