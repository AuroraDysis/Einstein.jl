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

        # Test 1D interior convolution with vector factor
        factors = [0.5, 0.6, 0.7, 0.8, 0.9]
        fdm_convolve_interior!(out_1d, in_1d, weights, factors)
        @test out_1d[2:4] ≈ [0.6, 0.7, 0.8]

        # Test 2D interior convolution with vector factor
        fdm_convolve_interior!(out_2d, in_2d, weights, factors)
        @test out_2d[2:4, :] ≈ [0.6 0.6; 0.7 0.7; 0.8 0.8]

        # Test Add mode for 1D interior convolution
        out_1d_add = ones(size(in_1d))
        fdm_convolve_interior!(out_1d_add, in_1d, weights, factor, ConvolveAadd())
        @test out_1d_add[2:4] ≈ [2.0, 2.0, 2.0]  # original ones + convolution result

        # Test Add mode for 2D interior convolution
        out_2d_add = ones(size(in_2d))
        fdm_convolve_interior!(out_2d_add, in_2d, weights, factor, ConvolveAadd())
        @test out_2d_add[2:4, :] ≈ [2.0 2.0; 2.0 2.0; 2.0 2.0]

        # Test Add mode for 1D interior convolution with vector factor
        out_1d_add = ones(size(in_1d))
        fdm_convolve_interior!(out_1d_add, in_1d, weights, factors, ConvolveAadd())
        @test out_1d_add[2:4] ≈ [1.6, 1.7, 1.8]

        # Test Add mode for 2D interior convolution with vector factor
        out_2d_add = ones(size(in_2d))
        fdm_convolve_interior!(out_2d_add, in_2d, weights, factors, ConvolveAadd())
        @test out_2d_add[2:4, :] ≈ [1.6 1.6; 1.7 1.7; 1.8 1.8]
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
        factors = [0.5, 0.6, 0.7, 0.8, 0.9]
        fdm_convolve_boundary!(out_1d, in_1d, left_weights, right_weights, factors)
        @test out_1d[1:2] ≈ [-1.5, 0.3]
        @test out_1d[4:5] ≈ [0.4, -5.4]

        # Test 2D boundary convolution
        in_2d = [1.0 2.0; 2.0 3.0; 3.0 4.0; 4.0 5.0; 5.0 6.0; 6.0 7.0]
        out_2d = similar(in_2d)
        left_weights = @SMatrix [1.0 -2.0; -0.5 0.5]
        right_weights = @SMatrix [-0.5 0.5; 1.0 -2.0]

        fdm_convolve_boundary!(out_2d, in_2d, left_weights, right_weights, factor)
        @test out_2d[1:2, :] ≈ [-3.0 -4.0; 0.5 0.5]
        @test out_2d[5:6, :] ≈ [0.5 0.5; -7.0 -8.0]

        # Test with vector factors
        fdm_convolve_boundary!(out_2d, in_2d, left_weights, right_weights, factors)
        @test out_2d[1:2, :] ≈ [-1.5 -2.0; 0.3 0.3]
        @test out_2d[5:6, :] ≈ [0.4 0.4; -6.3 -7.2]

        # Test Add mode for 1D boundary convolution
        out_1d_add = ones(size(in_1d))
        fdm_convolve_boundary!(out_1d_add, in_1d, left_weights, right_weights, factor, ConvolveAadd())
        @test out_1d_add[1:2] ≈ [-2.0, 1.5]  # original ones + convolution result
        @test out_1d_add[4:5] ≈ [1.5, -5.0]

        # Test Add mode with vector factors
        out_1d_add = ones(size(in_1d))
        fdm_convolve_boundary!(out_1d_add, in_1d, left_weights, right_weights, factors, ConvolveAadd())
        @test out_1d_add[1:2] ≈ [-0.5, 1.3]
        @test out_1d_add[4:5] ≈ [1.4, -4.4]

        # Test Add mode for 2D boundary convolution
        out_2d_add = ones(size(in_2d))
        fdm_convolve_boundary!(out_2d_add, in_2d, left_weights, right_weights, factor, ConvolveAadd())
        @test out_2d_add[1:2, :] ≈ [-2.0 -3.0; 1.5 1.5]
        @test out_2d_add[5:6, :] ≈ [1.5 1.5; -6.0 -7.0]

        # Test Add mode with vector factors for 2D
        out_2d_add = ones(size(in_2d))
        fdm_convolve_boundary!(out_2d_add, in_2d, left_weights, right_weights, factors, ConvolveAadd())
        @test out_2d_add[1:2, :] ≈ [-0.5 -1.0; 1.3 1.3]
        @test out_2d_add[5:6, :] ≈ [1.4 1.4; -5.3 -6.2]
    end
end
