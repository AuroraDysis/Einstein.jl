using TestItems
using Einstein.ChebyshevSuite, Test

@testitem "cheb1_integration_matrix" begin
    n = 10
    grid = ChebyshevGrid(n, 0.0, 1.0; kind=1)
    Q = cheb1_integration_matrix(n, 0.0, 1.0)
    @test size(Q) == (n, n)

    # Test scaling with interval change
    x_min, x_max = -2.0, 3.0
    Q_scaled = cheb1_integration_matrix(Float64, n, x_min, x_max)
    scale = (x_max - x_min) / 2
    Q_default = cheb1_integration_matrix(n)
    @test isapprox(Q_scaled, Q_default .* scale, rtol=1e-12)

    # Test integration on constant function: f = [1,1,...,1]
    f = ones(n)
    linear = Q * f
    for i in 1:n
        @test isapprox(linear[i], grid[i], rtol=1e-12)
    end
end
