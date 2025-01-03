@testset "cheb1_quadwts" begin
    # Test n=0 case
    @test cheb1_quadwts(0) == Float64[]

    # Test n=1 case
    @test cheb1_quadwts(1) ≈ [2.0]

    # Test n=5 case
    w5 = cheb1_quadwts(5)
    @test w5 ≈ [
        0.167781228466683,
        0.525552104866650,
        0.613333333333333,
        0.525552104866650,
        0.167781228466684,
    ]

    w6 = cheb1_quadwts(6)
    @test w6 ≈ [
        0.118661021381236,
        0.377777777777778,
        0.503561200840986,
        0.503561200840986,
        0.377777777777778,
        0.118661021381236,
    ]
end
