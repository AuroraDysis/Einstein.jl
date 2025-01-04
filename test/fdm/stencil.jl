@testitem "fdm_central" begin
    using GRSuite.FDMSuite, Test

    # Test 1st derivative, 2nd order accuracy
    @test fdm_central(1, 2) ≈ [-1//2, 0//1, 1//2]

    # Test 2nd derivative, 2nd order accuracy 
    @test fdm_central(2, 2) ≈ [1//1, -2//1, 1//1]

    # Test 1st derivative, 4th order accuracy
    @test fdm_central(1, 4) ≈ [1//12, -2//3, 0//1, 2//3, -1//12]

    # Test 6th derivative, 6th order accuracy
    @test fdm_central(6, 6) ≈ [13//240, -19//24, 87//16, -39//2, 323//8, -1023//20, 323//8, -39//2, 87//16, -19//24, 13//240]
end


@testitem "fdm_hermite" begin
    using GRSuite.FDMSuite, Test

    # Test 2nd derivative, 4th order accuracy
    D2, E2 = fdm_hermite(2, 4)
    @test D2 == [2//1, -4//1, 2//1]
    @test E2 == [1//2, 0//1, -1//2]

    # Test 2th derivative, 8nd order accuracy
    D4, E4 = fdm_hermite(2, 8)
    @test D4 == [7//54, 64//27, -5//1, 64//27, 7//54]
    @test E4 == [1//36, 8//9, 0//1, -8//9, -1//36]
end
