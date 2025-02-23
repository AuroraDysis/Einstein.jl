@testitem "fdm_dissipation_order" begin
    using Einstein.FiniteDifferenceSuite, Test
    # Test even accuracy orders
    @test fdm_dissipation_order(2) == 4
    @test fdm_dissipation_order(4) == 6
    @test fdm_dissipation_order(6) == 8

    # Test error for odd order
    @test_throws ArgumentError fdm_dissipation_order(3)
end

@testitem "fdm_dissipation_weights" begin
    using Einstein.FiniteDifferenceSuite, Test

    # Test 2nd order dissipation
    wts2 = fdm_dissipation_weights(2)
    wts2_expected = [1//4, -1//2, 1//4]
    @test wts2 == wts2_expected

    # Test 4th order dissipation
    wts4 = fdm_dissipation_weights(4)
    wts4_expected = [-1//16, 1//4, -3//8, 1//4, -1//16]
    @test wts4 == wts4_expected

    # Test 6th order dissipation
    wts6 = fdm_dissipation_weights(6)
    wts6_expected = [1//64, -3//32, 15//64, -5//16, 15//64, -3//32, 1//64]
    @test wts6 == wts6_expected
end
