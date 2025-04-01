using TestItems
using Einstein.ChebyshevSuite, Test

@testitem "integration_matrix" begin
    n = 10
    grid = ChebyshevGrid(n, 0.0, 1.0; kind=1)
    Q = integration_matrix(n, 0.0, 1.0)
    @test size(Q) == (n, n)

    # Test scaling with interval change
    lower_bound, upper_bound = -2.0, 3.0
    Q_scaled = integration_matrix(Float64, n, lower_bound, upper_bound)
    scale = (upper_bound - lower_bound) / 2
    Q_default = integration_matrix(n)
    @test isapprox(Q_scaled, Q_default .* scale, rtol=1e-12)

    # Test integration on constant function: f = [1,1,...,1]
    f = ones(n)
    linear = Q * f
    for i in 1:n
        @test isapprox(linear[i], grid[i], rtol=1e-12)
    end
end
