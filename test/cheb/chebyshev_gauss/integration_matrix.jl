@testitem "cheb_gauss_integration_matrix" begin
    n = 10
    tol = typetol(Float64)
    points = cheb_gauss_points(Float64, n, 0.0, 1.0)
    Q = cheb_gauss_integration_matrix(n, 1.0)
    @test size(Q) == (n, n)

    # Test scaling with interval change
    domain_width = 5.0
    Q_scaled = cheb_gauss_integration_matrix(Float64, n, domain_width)
    scale = domain_width / 2
    Q_default = cheb_gauss_integration_matrix(n)
    @test isapprox(Q_scaled, Q_default .* scale, atol=tol)

    # Test integration on constant function: f = [1,1,...,1]
    f = ones(n)
    linear = Q * f
    for i in 1:n
        @test isapprox(linear[i], points[i], atol=tol)
    end
end
