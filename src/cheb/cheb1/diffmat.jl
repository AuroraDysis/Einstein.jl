"""
    cheb1_diffmat(n::TI, k::TI=1) where {TI<:Integer}

Construct a Chebyshev differentiation matrix for points of the first kind.

# Description
Creates a matrix that maps function values at `n` Chebyshev points of the first kind 
to values of the `k`-th derivative of the interpolating polynomial at those points.

# Arguments
- `n::Integer`: Number of Chebyshev points
- `k::Integer=1`: Order of the derivative (default: 1)

# Returns
- Matrix that computes the `k`-th derivative when multiplied with function values

# Examples
```julia
# First derivative matrix for 5 points
D1 = cheb1_diffmat(5)

# Second derivative matrix for 5 points
D2 = cheb1_diffmat(5, 2)

# Apply to function values
vals = [1.0, 2.0, 3.0, 4.0, 5.0]
deriv_vals = D1 * vals
```

# References
- Based on the MATLAB Chebfun implementation
"""
function cheb1_diffmat(::Type{TR}, n::TI, k::TI=1) where {TR<:AbstractFloat,TI<:Integer}
    x = cheb1_pts(TR, n)               # First kind points.
    w = cheb1_barywts(TR, n)           # Barycentric weights.
    t = cheb1_angles(TR, n)            # acos(x).
    D = bary_diffmat(x, w, k, t)       # Construct matrix.
    return D
end

export cheb1_diffmat

@testset "cheb1_diffmat" begin
    @testset "n=5 first derivative" begin
        expected = [
            -4.97979656976556 7.20682929858878 -3.40260323340816 1.70130161670408 -0.525731112119134
            -1.05146222423827 -0.449027976579586 2.10292444847654 -0.850650808352040 0.248216560693358
            0.324919696232906 -1.37638192047117 -1.11022302462516e-16 1.37638192047117 -0.324919696232906
            -0.248216560693358 0.850650808352040 -2.10292444847654 0.449027976579586 1.05146222423827
            0.525731112119134 -1.70130161670408 3.40260323340816 -7.20682929858878 4.97979656976556
        ]
        result = cheb1_diffmat(Float64, 5)
        @test isapprox(result, expected, rtol=1e-12)
    end

    @testset "n=6 first derivative" begin
        expected = [
            -7.20976852010751 10.5558337350587 -5.27791686752937 3.04720672422855 -1.63299316185545 0.517638090205042
            -1.41421356237310 -0.707106781186547 3.04720672422855 -1.41421356237310 0.707106781186548 -0.218779599482357
            0.378937381963012 -1.63299316185545 -0.138700708242030 1.93185165257814 -0.757874763926024 0.218779599482357
            -0.218779599482357 0.757874763926024 -1.93185165257814 0.138700708242030 1.63299316185545 -0.378937381963012
            0.218779599482357 -0.707106781186548 1.41421356237310 -3.04720672422855 0.707106781186547 1.41421356237310
            -0.517638090205042 1.63299316185545 -3.04720672422855 5.27791686752937 -10.5558337350587 7.20976852010751
        ]
        result = cheb1_diffmat(Float64, 6)
        @test isapprox(result, expected, rtol=1e-12)
    end

    @testset "Type tests" begin
        @test eltype(cheb1_diffmat(Float32, 5)) == Float32
        @test eltype(cheb1_diffmat(Float64, 5)) == Float64
    end
end
