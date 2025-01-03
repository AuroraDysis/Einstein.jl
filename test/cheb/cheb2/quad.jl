using TestItems

@testitem "cheb2_quadwts" begin
    using PDESuite.ChebSuite, Test

    # Test n=0 case
    @test cheb2_quadwts(0) == Float64[]

    # Test n=1 case
    @test cheb2_quadwts(1) ≈ [2.0]

    # Test n=5 case
    w5 = cheb2_quadwts(5)
    @test w5 ≈ [
        0.0666666666666667,
        0.533333333333333,
        0.800000000000000,
        0.533333333333333,
        0.0666666666666667,
    ]

    w6 = cheb2_quadwts(6)
    @test w6 ≈ [
        0.0400000000000000,
        0.360743041200011,
        0.599256958799989,
        0.599256958799989,
        0.360743041200011,
        0.0400000000000000,
    ]
end
