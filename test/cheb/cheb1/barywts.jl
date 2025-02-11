using TestItems

@testitem "cheb1_barycentric_weights" begin
    using Einstein.ChebyshevSuite, Test

    # Test n=0 case
    @test cheb1_barycentric_weights(0) == Float64[]

    # Test n=1 case
    @test cheb1_barycentric_weights(1) ≈ [1.0]

    # Test n=5 case
    w5 = cheb1_barycentric_weights(5)
    @test w5 ≈ [
        0.309016994374947, -0.809016994374948, 1.0, -0.809016994374948, 0.309016994374947
    ]

    w6 = cheb1_barycentric_weights(6)
    @test w6 ≈ [
        -0.258819045102521,
        0.707106781186548,
        -0.965925826289068,
        0.965925826289068,
        -0.707106781186548,
        0.258819045102521,
    ]
end
