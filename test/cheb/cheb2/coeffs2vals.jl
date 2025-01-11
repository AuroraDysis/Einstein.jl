using TestItems

@testitem "cheb2_coeffs2vals" begin
    using LinearAlgebra, Einstein.ChebSuite, Test

    # Set tolerance
    tol = 100 * eps()

    # Test single coefficient conversion
    c = [sqrt(2)]
    v = cheb2_coeffs2vals(c)
    @test c ≈ v

    c = collect(1.0:5.0)
    vTrue = [3; -4 + sqrt(2); 3; -4 - sqrt(2); 15]
    v = cheb2_coeffs2vals(c)
    @test norm(v - vTrue, Inf) < tol

    c = collect(1.0:6.0)
    vTrue = [-3; 7 / 2; -(11 / 2) + sqrt(5); 7 / 2; -(11 / 2) - sqrt(5); 21]
    v = cheb2_coeffs2vals(c)
    @test norm(v - vTrue, Inf) < tol

    # Test symmetry preservation
    c = kron(ones(10), Matrix{Float64}(I, 2, 2))
    c1 = @view c[:, 1]
    c2 = @view c[:, 2]
    v1 = cheb2_coeffs2vals(c1)
    v2 = cheb2_coeffs2vals(c2)
    @test norm(v1 - reverse(v1), Inf) ≈ 0
    @test norm(v2 + reverse(v2), Inf) ≈ 0
end
