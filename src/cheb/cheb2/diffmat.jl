"""
    cheb2_diffmat([TR=Float64], n::TI, k::TI=1) where {TR<:AbstractFloat,TI<:Integer}

Construct a Chebyshev differentiation that maps function values at `n` Chebyshev points of the 2nd kind 
to values of the `k`-th derivative of the interpolating polynomial at those points.

# Arguments
- `TR`: Element type (defaults to Float64)
- `n::Integer`: Number of Chebyshev points
- `k::Integer=1`: Order of the derivative (default: 1)

# References
- [chebfun/@chebcolloc2/chebcolloc2.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc2/chebcolloc2.m)
"""
function cheb2_diffmat(::Type{TR}, n::TI, k::TI=1) where {TR<:AbstractFloat,TI<:Integer}
    x = cheb2_pts(TR, n)               # First kind points.
    w = cheb2_barywts(TR, n)           # Barycentric weights.
    t = cheb2_angles(TR, n)            # acos(x).
    D = bary_diffmat(x, w, k, t)       # Construct matrix.
    return D
end

function cheb2_diffmat(n::TI, k::TI=1) where {TI<:Integer}
    return cheb2_diffmat(Float64, n, k)
end

@testset "cheb2_diffmat" begin
    tol = 100 * eps()
    # Test case for n=5
    D5 = cheb2_diffmat(Float64, 5)
    D5_expected = [
        -5.5 6.82842712474619 -2.0 1.17157287525381 -0.5
        -1.70710678118655 0.707106781186548 1.41421356237310 -0.707106781186548 0.292893218813453
        0.5 -1.41421356237310 0.0 1.41421356237310 -0.5
        -0.292893218813453 0.707106781186548 -1.41421356237310 -0.707106781186548 1.70710678118655
        0.5 -1.17157287525381 2.0 -6.82842712474619 5.5
    ]
    @test isapprox(D5, D5_expected, rtol=tol)

    # Test case for n=6
    D6 = cheb2_diffmat(Float64, 6)
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
