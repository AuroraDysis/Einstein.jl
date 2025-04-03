@testitem "cheb_gauss_vals2coeffs_matrix, cheb_gauss_coeffs2vals_matrix" begin
    for type in [Float64, BigFloat]
        tol = typetol(type)

        @testset "matrix properties" begin
            for n in 1:10
                A = cheb_gauss_vals2coeffs_matrix(type, n)
                S = cheb_gauss_coeffs2vals_matrix(type, n)

                # Test matrix dimensions
                @test size(A) == (n, n)
                @test size(S) == (n, n)

                # Test that matrices are inverses of each other
                @test isapprox(A * S, Matrix(I, n, n), atol=tol)
                @test isapprox(S * A, Matrix(I, n, n), atol=tol)

                # Test that matrices are symmetric
                @test isapprox(A, A', atol=tol)
                @test isapprox(S, S', atol=tol)
            end
        end

        @testset "polynomial interpolation" begin
            for n in 1:10
                # Create a polynomial function
                f(x) = 3x^2 + 2x - 1
                
                # Get points and evaluate function
                x = cheb_gauss_points(type, n)
                f_vals = f.(x)
                
                # Transform to coefficients and back
                A = cheb_gauss_vals2coeffs_matrix(type, n)
                S = cheb_gauss_coeffs2vals_matrix(type, n)
                
                f_coeffs = A * f_vals
                f_vals_reconstructed = S * f_coeffs
                
                # Test reconstruction
                @test isapprox(f_vals, f_vals_reconstructed, atol=tol)
            end
        end

        @testset "trigonometric interpolation" begin
            for n in 1:10
                # Create a trigonometric function
                f(x) = sin(π * x) * cos(2π * x)
                
                # Get points and evaluate function
                x = cheb_gauss_points(type, n)
                f_vals = f.(x)
                
                # Transform to coefficients and back
                A = cheb_gauss_vals2coeffs_matrix(type, n)
                S = cheb_gauss_coeffs2vals_matrix(type, n)
                
                f_coeffs = A * f_vals
                f_vals_reconstructed = S * f_coeffs
                
                # Test reconstruction
                @test isapprox(f_vals, f_vals_reconstructed, atol=tol)
            end
        end
    end
end