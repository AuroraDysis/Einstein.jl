using TestItems

@testitem "cheb1_barywts" begin
    using Einstein.ChebSuite, Test

    # Test n=0 case
    @test cheb1_barywts(0) == Float64[]

    # Test n=1 case
    @test cheb1_barywts(1) ≈ [1.0]

    # Test n=5 case
    w5 = cheb1_barywts(5)
    @test w5 ≈ [
        0.309016994374947, -0.809016994374948, 1.0, -0.809016994374948, 0.309016994374947
    ]

    w6 = cheb1_barywts(6)
    @test w6 ≈ [
        -0.258819045102521,
        0.707106781186548,
        -0.965925826289068,
        0.965925826289068,
        -0.707106781186548,
        0.258819045102521,
    ]
end
