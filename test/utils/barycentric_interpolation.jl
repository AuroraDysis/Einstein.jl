@testitem "Barycentric Interpolation" begin
    using LinearAlgebra

    tol = typetol(Float64)

    @testset "Barycentric Weights" begin
        # Test equispaced nodes
        @testset "Equispaced nodes" begin
            # Test using explicit array
            x_eq = collect(range(-1.0, 1.0, 5))
            w_eq1 = barycentric_weights(x_eq)
            
            # Test using range
            x_range = range(-1.0, 1.0, 5)
            w_eq2 = barycentric_weights(x_range)
            
            # Test using order specification
            w_eq3 = barycentric_weights(Float64, 4)  # order = length(x) - 1 = 4
            
            # All should be the same
            @test w_eq1 ≈ w_eq2
            @test w_eq1 ≈ w_eq3
            @test norm(w_eq1, Inf) ≈ 1.0  # Weights are normalized
            
            # Check alternating signs pattern for equispaced nodes
            for i in 1:length(w_eq1)-1
                @test sign(w_eq1[i]) == -sign(w_eq1[i+1])
            end
        end
        
        # Test non-equispaced nodes (Chebyshev points)
        @testset "Non-equispaced nodes (Chebyshev)" begin
            n = 5
            # Chebyshev points of the second kind
            cheb_points = cheb_lobatto_points(n)
            w_cheb = barycentric_weights(cheb_points)
            
            # Check normalization
            @test norm(w_cheb, Inf) ≈ 1.0
            
            # Check that weights are correct for Chebyshev points
            # For Chebyshev points, weights should have alternating signs with
            # special values at endpoints
            for i in 1:length(w_cheb)-1
                if i == 1 || i == length(w_cheb)
                    @test abs(w_cheb[i]) ≈ 0.5 atol=tol
                end
                @test sign(w_cheb[i]) == -sign(w_cheb[i+1])
            end
        end
        
        # Test for errors with invalid inputs
        @testset "Error handling" begin
            # Single point
            @test_throws ArgumentError barycentric_weights([1.0])
            
            # Zero or negative order
            @test_throws ArgumentError barycentric_weights(Float64, 0)
            @test_throws ArgumentError barycentric_weights(Float64, -1)
            
            # Points not distinct (all the same)
            @test_throws ArgumentError barycentric_weights([1.0, 1.0, 1.0])
        end
    end

    @testset "Barycentric Interpolation" begin
        @testset "Polynomial exact recovery" begin
            # Test polynomial recovery with equispaced nodes
            @testset "Equispaced nodes" begin
                # Define a polynomial: p(x) = 1 + 2x + 3x² + 4x³
                f(x) = 1.0 + 2.0*x + 3.0*x^2 + 4.0*x^3
                
                # Interpolation points
                x = range(-1.0, 1.0, 8)
                y = f.(x)
                
                # Create weights and interpolator
                weights = barycentric_weights(x)
                
                # Test at interpolation points (should be exact)
                for (xi, yi) in zip(x, y)
                    @test barycentric_interpolate(xi, collect(x), y, weights) ≈ yi
                end
                
                # Test at non-interpolation points (should be very close for polynomial)
                test_points = range(-0.95, 0.95, 10)
                for t in test_points
                    interp_val = barycentric_interpolate(t, collect(x), y, weights)
                    exact_val = f(t)
                    @test interp_val ≈ exact_val atol=tol
                end
                
                # Test using BarycentricInterpolation struct
                itp = BarycentricInterpolation(collect(x), weights)
                for t in test_points
                    @test itp(y, t) ≈ f(t) atol=tol
                end
            end
            
            # Test polynomial recovery with non-equispaced nodes (Chebyshev)
            @testset "Non-equispaced nodes (Chebyshev)" begin
                # Define a polynomial: p(x) = 1 + 2x + 3x² + 4x³ + 5x⁴
                f(x) = 1.0 + 2.0*x + 3.0*x^2 + 4.0*x^3 + 5.0*x^4
                
                # Chebyshev points (more efficient for high-degree polynomials)
                n = 9
                cheb_points = [cos((j * π) / (n - 1)) for j in 0:n-1]
                sort!(cheb_points)  # Ensure sorted order
                
                y = f.(cheb_points)
                
                # Create weights and interpolator
                weights = barycentric_weights(cheb_points)
                
                # Test at interpolation points (should be exact)
                for (xi, yi) in zip(cheb_points, y)
                    @test barycentric_interpolate(xi, cheb_points, y, weights) ≈ yi
                end
                
                # Test at non-interpolation points
                test_points = range(-0.95, 0.95, 10)
                for t in test_points
                    interp_val = barycentric_interpolate(t, cheb_points, y, weights)
                    exact_val = f(t)
                    @test interp_val ≈ exact_val atol=tol
                end
                
                # Test BarycentricInterpolation struct
                itp = BarycentricInterpolation(cheb_points, weights)
                for t in test_points
                    @test itp(y, t) ≈ f(t) atol=tol
                end
            end
        end
        
        @testset "Trigonometric function interpolation" begin
            # Test equispaced vs non-equispaced nodes for interpolating sin(5x)
            # This will demonstrate the superiority of non-equispaced nodes for 
            # oscillatory functions
            f(x) = sin(5.0 * x)
            
            # Compare equispaced and Chebyshev points for different numbers of points
            for n in [11, 21, 31]
                # Equispaced points
                x_eq = range(-1.0, 1.0, n)
                y_eq = f.(x_eq)
                w_eq = barycentric_weights(x_eq)
                
                # Chebyshev points
                cheb_points = [cos((j * π) / (n - 1)) for j in 0:n-1]
                sort!(cheb_points)  # Ensure sorted order
                y_cheb = f.(cheb_points)
                w_cheb = barycentric_weights(cheb_points)
                
                # Test points (dense grid for good error estimation)
                test_points = range(-0.95, 0.95, 100)
                
                # Maximum error for equispaced
                max_err_eq = maximum(abs.(f.(test_points) .- 
                    [barycentric_interpolate(t, collect(x_eq), y_eq, w_eq) for t in test_points]))
                
                # Maximum error for Chebyshev 
                max_err_cheb = maximum(abs.(f.(test_points) .- 
                    [barycentric_interpolate(t, cheb_points, y_cheb, w_cheb) for t in test_points]))
                
                # Chebyshev points should generally give better results for this function
                # Especially as n increases
                if n > 20
                    @test max_err_cheb < max_err_eq
                end
            end
        end
        
        @testset "Error handling" begin
            x = [-1.0, 0.0, 1.0]
            y = [1.0, 0.0, 1.0]
            w = barycentric_weights(x)
            
            # Test bounds checking in BarycentricInterpolation
            itp = BarycentricInterpolation(x, w)
            @test_throws ArgumentError itp(y, -1.5)  # out of range
            @test_throws ArgumentError itp(y, 1.5)   # out of range
            
            # Points, values and weights size mismatch
            @test_throws ArgumentError barycentric_interpolate(0.5, x, [1.0, 0.0], w)
            @test_throws ArgumentError barycentric_interpolate(0.5, x, y, [1.0, 2.0])
            
            # Points too few
            @test_throws ArgumentError barycentric_interpolate(0.5, [1.0], [1.0], [1.0])
        end
    end
end
