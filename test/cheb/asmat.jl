@testitem "cheb_vals2coeffs_matrix, cheb_coeffs2vals_matrix" begin
    using LinearAlgebra

    tol = typetol(Float64)
    grid_list = [ChebyshevGaussGrid, ChebyshevLobattoGrid]

    for grid in grid_list
        @testset "$grid, matrix properties" begin
            for n in 1:10
                A = grid.vals2coeffs_matrix(n)
                S = grid.coeffs2vals_matrix(n)

                # Test matrix dimensions
                @test size(A) == (n, n)
                @test size(S) == (n, n)

                # Test that matrices are inverses of each other
                @test isapprox(A * S, Matrix(I, n, n), atol=tol)
                @test isapprox(S * A, Matrix(I, n, n), atol=tol)
            end
        end

        @testset "$grid, polynomial interpolation" begin
            for n in 1:10
                # Create a polynomial function
                f(x) = 3x^2 + 2x - 1

                # Get points and evaluate function
                x = grid.points(n)
                f_vals = f.(x)

                # Transform to coefficients and back
                A = grid.vals2coeffs_matrix(n)
                S = grid.coeffs2vals_matrix(n)

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
                x = grid.points(n)
                f_vals = f.(x)

                # Transform to coefficients and back
                A = grid.vals2coeffs_matrix(n)
                S = grid.coeffs2vals_matrix(n)

                f_coeffs = A * f_vals
                f_vals_reconstructed = S * f_coeffs

                # Test reconstruction
                @test isapprox(f_vals, f_vals_reconstructed, atol=tol)
            end
        end
    end
end
