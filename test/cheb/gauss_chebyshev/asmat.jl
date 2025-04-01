using TestItems

@testitem "GaussChebyshev - vals2coeffs_matrix, coeffs2vals_matrix" begin
    using ApproxFun

    for type in [Float64, BigFloat]
        tol = typetol(type)

        @testset "polynomial" begin
            dom = -one(type) .. one(type)
            f = Fun(x -> 3x^2 + 2x - 1, Chebyshev(dom))
            f_coeffs = coefficients(f)
            n = ncoefficients(f)

            x = GaussChebyshev.points(type, n)
            f_vals = f.(x)
            A = GaussChebyshev.vals2coeffs_matrix(type, n)
            S = GaussChebyshev.coeffs2vals_matrix(type, n)

            @test isapprox(A * f_vals, f_coeffs, atol=tol)
            @test isapprox(S * f_coeffs, f_vals, atol=tol)
        end

        @testset "trigonometric" begin
            dom = -one(type) .. one(type)
            f = Fun(x -> sin(π * x) * cos(2π * x), Chebyshev(dom))
            f_coeffs = coefficients(f)
            n = ncoefficients(f)

            x = GaussChebyshev.points(type, n)
            f_vals = f.(x)
            A = GaussChebyshev.vals2coeffs_matrix(type, n)
            S = GaussChebyshev.coeffs2vals_matrix(type, n)

            @test isapprox(A * f_vals, f_coeffs, atol=tol)
            @test isapprox(S * f_coeffs, f_vals, atol=tol)
        end
    end
end