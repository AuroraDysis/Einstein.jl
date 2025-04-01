using TestItems

@testitem "GaussChebyshevGrid - differentiation_matrix" begin
    # n = 5
    expected5 = [
        -4.97979656976556 7.20682929858878 -3.40260323340816 1.70130161670408 -0.525731112119134
        -1.05146222423827 -0.449027976579586 2.10292444847654 -0.850650808352040 0.248216560693358
        0.324919696232906 -1.37638192047117 -1.11022302462516e-16 1.37638192047117 -0.324919696232906
        -0.248216560693358 0.850650808352040 -2.10292444847654 0.449027976579586 1.05146222423827
        0.525731112119134 -1.70130161670408 3.40260323340816 -7.20682929858878 4.97979656976556
    ]
    result5 = GaussChebyshevGrid.differentiation_matrix(Float64, 5)
    @test isapprox(result5, expected5, rtol=1e-12)

    # n = 6
    expected6 = [
        -7.20976852010751 10.5558337350587 -5.27791686752937 3.04720672422855 -1.63299316185545 0.517638090205042
        -1.41421356237310 -0.707106781186547 3.04720672422855 -1.41421356237310 0.707106781186548 -0.218779599482357
        0.378937381963012 -1.63299316185545 -0.138700708242030 1.93185165257814 -0.757874763926024 0.218779599482357
        -0.218779599482357 0.757874763926024 -1.93185165257814 0.138700708242030 1.63299316185545 -0.378937381963012
        0.218779599482357 -0.707106781186548 1.41421356237310 -3.04720672422855 0.707106781186547 1.41421356237310
        -0.517638090205042 1.63299316185545 -3.04720672422855 5.27791686752937 -10.5558337350587 7.20976852010751
    ]
    result6 = GaussChebyshevGrid.differentiation_matrix(Float64, 6)
    @test isapprox(result6, expected6, rtol=1e-12)

    # Type tests
    @test eltype(GaussChebyshevGrid.differentiation_matrix(Float32, 5)) == Float32
    @test eltype(GaussChebyshevGrid.differentiation_matrix(Float64, 5)) == Float64
end
