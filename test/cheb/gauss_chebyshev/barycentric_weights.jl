using TestItems

@testitem "cheb_gauss_barycentric_weights" begin
    weights_0 = Float64[]
    weights_1 = [1.0]
    weights_5 = [
        0.309016994374947, -0.809016994374948, 1.0, -0.809016994374948, 0.309016994374947
    ]
    weights_6 = [
        -0.258819045102521,
        0.707106781186548,
        -0.965925826289068,
        0.965925826289068,
        -0.707106781186548,
        0.258819045102521,
    ]

    @test cheb_gauss_barycentric_weights(0) ≈ weights_0
    @test cheb_gauss_barycentric_weights(1) ≈ weights_1
    @test cheb_gauss_barycentric_weights(5) ≈ weights_5
    @test cheb_gauss_barycentric_weights(6) ≈ weights_6
end
