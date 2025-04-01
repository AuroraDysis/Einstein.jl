using TestItems
using Einstein.ChebyshevSuite, Test

@testitem "GaussChebyshevGrid - integration_matrix" begin
    n = 10
    points = GaussChebyshevGrid.points(Float64, n, 0.0, 1.0)
    Q = GaussChebyshevGrid.integration_matrix(n, 0.0, 1.0)
    @test size(Q) == (n, n)

    # Test scaling with interval change
    lower_bound, upper_bound = -2.0, 3.0
    Q_scaled = GaussChebyshevGrid.integration_matrix(Float64, n, lower_bound, upper_bound)
    scale = (upper_bound - lower_bound) / 2
    Q_default = GaussChebyshevGrid.integration_matrix(n)
    @test isapprox(Q_scaled, Q_default .* scale, rtol=1e-12)

    # Test integration on constant function: f = [1,1,...,1]
    f = ones(n)
    linear = Q * f
    for i in 1:n
        @test isapprox(linear[i], points[i], rtol=1e-12)
    end
end
