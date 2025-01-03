using TestItems

@testitem "cheb2_barywts" begin
    using PDESuite.ChebSuite, Test

    # Test n=1 case
    w1 = cheb2_barywts(Float64, 1)
    @test w1 ≈ [1.0]

    # Test n=2 case
    w2 = cheb2_barywts(Float64, 2)
    @test w2 ≈ [-0.500000000000000, 0.500000000000000]

    # Test n=5 case
    w5 = cheb2_barywts(Float64, 5)
    @test w5 ≈ [0.500000000000000, -1.0, 1.0, -1.0, 0.500000000000000]

    w6 = cheb2_barywts(Float64, 6)
    @test w6 ≈ [-0.500000000000000, 1.0, -1.0, 1.0, -1.0, 0.500000000000000]
end
