using TestItems

@testitem "cheb_cumsum" begin
    using LinearAlgebra, PDESuite.ChebSuite, Test

    tol = 100 * eps()

    @testset "cheb1" begin
        n = 15
        f = cos.(cheb1_pts(n))
        f_coeffs = cheb1_vals2coeffs(f)
        If_coeffs = cheb_cumsum(f_coeffs)
        If = sin.(cheb1_pts(n)) .- sin(-1) # sin(x) - sin(-1) is the antiderivative of cos(x)
        If_coeffs_true = cheb1_vals2coeffs(If)
        @test norm(If_coeffs[1:(end - 1)] .- If_coeffs_true, Inf) < tol
    end

    @testset "cheb2" begin
        n = 15
        f = cos.(cheb2_pts(n))
        f_coeffs = cheb2_vals2coeffs(f)
        If_coeffs = cheb_cumsum(f_coeffs)
        If = sin.(cheb2_pts(n)) .- sin(-1) # sin(x) - sin(-1) is the antiderivative of cos(x)
        If_coeffs_true = cheb2_vals2coeffs(If)
        @test norm(If_coeffs[1:(end - 1)] .- If_coeffs_true, Inf) < tol
    end
end
