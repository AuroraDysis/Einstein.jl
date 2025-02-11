using TestItems

@testitem "cheb1_analysis_matrix, cheb1_synthesis_matrix" begin
    using ApproxFun, Einstein.ChebyshevSuite, Test

    for type in [Float64, BigFloat]
        tol = 1000 * eps(type)

        @testset "polynomial" begin
            dom = -one(type) .. one(type)
            f = Fun(x -> 3x^2 + 2x - 1, Chebyshev(dom))
            f_coeffs = coefficients(f)
            n = ncoefficients(f)

            x = cheb1_points(type, n)
            f_vals = f.(x)
            A = cheb1_analysis_matrix(type, n)
            S = cheb1_synthesis_matrix(type, n)

            @test isapprox(A * f_vals, f_coeffs, atol=tol)
            @test isapprox(S * f_coeffs, f_vals, atol=tol)
        end

        @testset "trigonometric" begin
            dom = -one(type) .. one(type)
            f = Fun(x -> sin(π * x) * cos(2π * x), Chebyshev(dom))
            f_coeffs = coefficients(f)
            n = ncoefficients(f)

            x = cheb1_points(type, n)
            f_vals = f.(x)
            A = cheb1_analysis_matrix(type, n)
            S = cheb1_synthesis_matrix(type, n)

            @test isapprox(A * f_vals, f_coeffs, atol=tol)
            @test isapprox(S * f_coeffs, f_vals, atol=tol)
        end
    end
end