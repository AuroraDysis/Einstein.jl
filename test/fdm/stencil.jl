@testitem "fdm_central" begin
    using Einstein.FDMSuite, Test

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
    stencil = [
        13//240,
        -19//24,
        87//16,
        -39//2,
        323//8,
        -1023//20,
        323//8,
        -39//2,
        87//16,
        -19//24,
        13//240,
    ]
    @test fdm_central(6, 6) ≈ stencil
    for type in [Float64, BigFloat]
        @test fdm_central(type, 6, 6) ≈ type.(stencil)
    end
end

@testitem "fdm_hermite" begin
    using Einstein.FDMSuite, Test

    # Test 2nd derivative, 4th order accuracy
    D2, E2 = fdm_hermite(2, 4)
    @test D2 == [2//1, -4//1, 2//1]
    @test E2 == [1//2, 0//1, -1//2]

    # Test 2th derivative, 8nd order accuracy
    D4, E4 = fdm_hermite(2, 8)
    @test D4 == [7//54, 64//27, -5//1, 64//27, 7//54]
    @test E4 == [1//36, 8//9, 0//1, -8//9, -1//36]
end

@testitem "fdm_extrapwts" begin
    using Einstein.FDMSuite, Test

    @test fdm_extrapwts_left(4) == [4, -6, 4, -1]
    @test fdm_extrapwts_right(4) == [-1, 4, -6, 4]

    @test fdm_extrapwts_left(5) == [5, -10, 10, -5, 1]
    @test fdm_extrapwts_right(5) == [1, -5, 10, -10, 5]
end

@testitem "fdm_boundwts" begin
    using Einstein.FDMSuite, Test

    Dl14, Dr14 = fdm_boundwts(1, 4)
    @test Dl14[:, 1] == [-25//12, 4//1, -3//1, 4//3, -1//4]
    @test Dl14[:, 2] == [-1//4, -5//6, 3//2, -1//2, 1//12]
    @test Dr14[:, 1] == [-1//12, 1//2, -3//2, 5//6, 1//4]
    @test Dr14[:, 2] == [1//4, -4//3, 3//1, -4//1, 25//12]

    Dl24, Dr24 = fdm_boundwts(2, 4)
    @test Dl24[:, 1] == [15//4, -77//6, 107//6, -13//1, 61//12, -5//6]
    @test Dl24[:, 2] == [5//6, -5//4, -1//3, 7//6, -1//2, 1//12]
    @test Dr24[:, 1] == [1//12, -1//2, 7//6, -1//3, -5//4, 5//6]
    @test Dr24[:, 2] == [-5//6, 61//12, -13//1, 107//6, -77//6, 15//4]
end
