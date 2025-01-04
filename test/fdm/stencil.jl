@testitem "fdm_central" begin
    using GRSuite.FDMSuite, Test

    # Test 1st derivative, 2nd order accuracy
    stencil = [-1//2, 0//1, 1//2]
    @test fdm_central(1, 2) == stencil
    for type in [Float64, BigFloat]
        @test fdm_central(type, 1, 2) ≈ type.(stencil)
    end

    # Test 2nd derivative, 2nd order accuracy 
    stencil = [1//1, -2//1, 1//1]
    @test fdm_central(2, 2) ≈ stencil
    for type in [Float64, BigFloat]
        @test fdm_central(type, 2, 2) ≈ type.(stencil)
    end

    # Test 1st derivative, 4th order accuracy
    stencil = [1//12, -2//3, 0//1, 2//3, -1//12]
    @test fdm_central(1, 4) == stencil
    for type in [Float64, BigFloat]
        @test fdm_central(type, 1, 4) ≈ type.(stencil)
    end

    # Test 6th derivative, 6th order accuracy
    stencil = [13//240, -19//24, 87//16, -39//2, 323//8, -1023//20, 323//8, -39//2, 87//16, -19//24, 13//240]
    @test fdm_central(6, 6) ≈ stencil
    for type in [Float64, BigFloat]
        @test fdm_central(type, 6, 6) ≈ type.(stencil)
    end
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
