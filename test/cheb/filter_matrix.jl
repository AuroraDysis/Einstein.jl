@testitem "Chebyshev Filter Matrix" begin
    using LinearAlgebra

    n = 50
    α = 36
    p = 32
    tol = typetol(Float64)

    grid_list = [ChebyshevGaussGrid, ChebyshevLobattoGrid]
    negative_sum_trick_list = [true, false]
    for grid in grid_list
        for negative_sum_trick in negative_sum_trick_list
            @testset "Negative sum trick: $negative_sum_trick" begin
                x = grid.points(n)
                weights = cheb_series_filter_weights_exp(n, α, p)
                S = grid.coeffs2vals_matrix(n)
                A = grid.vals2coeffs_matrix(n)
                F = cheb_filter_matrix(weights, S, A; negative_sum_trick=negative_sum_trick)

                # Check dimensions
                @test size(F) == (n, n)

                # Check filter on constant function, should be identity
                f_constant(x) = 1.0
                f_constant_vals = f_constant.(x)
                f_constant_filtered = F * f_constant_vals
                @test isapprox(f_constant_filtered, f_constant_vals, atol=tol)

                # Check filter on smooth function
                f_vals = sin.(π * x)
                f_filtered = F * f_vals
                @test isapprox(f_filtered, f_vals, atol=tol)
            end
        end
    end
end
