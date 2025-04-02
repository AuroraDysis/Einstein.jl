using TestItems

@testitem "cheb_rect_integration_matrix" begin
    n = 4
    tol = typetol(Float64)
    intmat = cheb_rect_integration_matrix(n)
    intmat_ana = [
        -1.73472347597681e-17 -2.08166817117217e-17 2.08166817117217e-17 3.46944695195361e-18
        0.204861111111111 0.326388888888889 -0.0486111111111111 0.0173611111111111
        0.0937500000000000 0.937500000000000 0.562500000000000 -0.0937500000000000
        0.111111111111111 0.888888888888889 0.888888888888889 0.111111111111111
    ]
    @test isapprox(intmat, intmat_ana; atol=tol)

    n = 5
    intmat = cheb_rect_integration_matrix(n)
    intmat_ana = [
        -3.46944695195361e-18 1.38777878078145e-17 2.08166817117217e-17 -6.93889390390723e-18 0
        0.119403559372885 0.190063432038124 -0.0242640687119285 0.0132867367414871 -0.00559644062711508
        0.0333333333333333 0.620220057259940 0.400000000000000 -0.0868867239266071 0.0333333333333333
        0.0722631072937817 0.520046596591846 0.824264068711929 0.343269901295209 -0.0527368927062182
        0.0666666666666667 0.533333333333333 0.800000000000000 0.533333333333333 0.0666666666666667
    ]

    @test isapprox(intmat, intmat_ana; atol=tol)

    @testset "Compare direct vs coefficient-based integration" begin
        for n in [4, 8, 16, 32]
            # Test on standard domain [-1,1]
            I1 = cheb_rect_integration_matrix(n)
            I2 = cheb_lobatto_integration_matrix(Float64, n)
            @test isapprox(I1, I2; atol=tol)

            # Test on mapped domain [0,π]
            I1 = cheb_rect_integration_matrix(n, 0.0, Float64(π))
            I2 = cheb_lobatto_integration_matrix(Float64, n, 0.0, Float64(π))
            @test isapprox(I1, I2; atol=tol)
        end
    end
end

@testitem "cheb_rect_integration_matrix - analytical" begin
    @testset "Standard domain [-1,1]" begin
        n = 32
        x = cheb_lobatto_points(Float64, n)
        intmat = cheb_rect_integration_matrix(n)

        # Test 1: Polynomial integration
        f = @. x^3  # f(x) = x³
        F_numeric = intmat * f   # Should give (x⁴ - 1) / 4
        F_exact = @. (x^4 - 1) / 4
        @test isapprox(F_numeric, F_exact; atol=tol)

        # Test 2: Trigonometric integration
        f = @. sin(π * x)
        F_numeric = intmat * f   # Should give -(cos(πx)+1)/π
        F_exact = @. -(cos(π * x) + 1) / π
        @test isapprox(F_numeric, F_exact; atol=tol)
    end

    @testset "Mapped domain [0,π]" begin
        n = 32
        intmat = cheb_rect_integration_matrix(n, 0.0, Float64(π))
        x = GaussChebyshevLobattoGrid.points(n, 0.0, Float64(π))

        # Test: Integration of sin(x) from 0 to x
        f = sin.(x)
        F_numeric = intmat * f   # Should give -cos(x) + 1
        F_exact = @. -cos(x) + 1
        @test isapprox(F_numeric, F_exact; atol=tol)

        # Test: Integration of x*cos(x)
        f = @. x * cos(x)
        F_numeric = intmat * f   # Should give x*sin(x) - sin(x)
        F_exact = @. x * sin(x) + cos(x) - 1
        @test isapprox(F_numeric, F_exact; atol=tol)
    end
end
