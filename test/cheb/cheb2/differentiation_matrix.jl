using TestItems

@testitem "differentiation_matrix" begin
    using Einstein.ChebyshevSuite, Test

    tol = 100 * eps()
    # Test case for n=5
    D5 = differentiation_matrix(Float64, 5)
    D5_expected = [
        -5.5 6.82842712474619 -2.0 1.17157287525381 -0.5
        -1.70710678118655 0.707106781186548 1.41421356237310 -0.707106781186548 0.292893218813453
        0.5 -1.41421356237310 0.0 1.41421356237310 -0.5
        -0.292893218813453 0.707106781186548 -1.41421356237310 -0.707106781186548 1.70710678118655
        0.5 -1.17157287525381 2.0 -6.82842712474619 5.5
    ]
    @test isapprox(D5, D5_expected, rtol=tol)

    # Test case for n=6
    D6 = differentiation_matrix(Float64, 6)
    D6_expected = [
        -8.5 10.4721359549996 -2.89442719099992 1.52786404500042 -1.10557280900008 0.5
        -2.61803398874990 1.17082039324994 2.0 -0.894427190999916 0.618033988749895 -0.276393202250021
        0.723606797749979 -2.0 0.170820393249937 1.61803398874990 -0.894427190999916 0.381966011250105
        -0.381966011250105 0.894427190999916 -1.61803398874990 -0.170820393249937 2.0 -0.723606797749979
        0.276393202250021 -0.618033988749895 0.894427190999916 -2.0 -1.17082039324994 2.61803398874990
        -0.5 1.10557280900008 -1.52786404500042 2.89442719099992 -10.4721359549996 8.5
    ]
    @test isapprox(D6, D6_expected, rtol=tol)
end
