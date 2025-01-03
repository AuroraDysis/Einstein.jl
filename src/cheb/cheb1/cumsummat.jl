"""
    cheb1_cumsummat([TR=Float64], n::Integer) where {TR<:AbstractFloat}
    cheb1_cumsummat([TR=Float64], n::Integer, x_min::TR, x_max::TR) where {TR<:AbstractFloat}

Compute Chebyshev integration matrix that maps function values
at `n` Chebyshev points of the 1st kind to values of the integral of the interpolating
polynomial at those points, with the convention that the first value is zero.

# References
- [chebfun/@chebcolloc1/chebcolloc1.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc1/chebcolloc1.m)
"""
function cheb1_cumsummat(::Type{TR}, n::Integer) where {TR<:AbstractFloat}
    A = cheb1_amat(TR, n)
    S = cheb1_smat(TR, n)
    B = cheb_coeffs_cumsummat(TR, n)
    Q = S * B * A
    return Q
end

function cheb1_cumsummat(n::Integer)
    return cheb1_cumsummat(Float64, n)
end

function cheb1_cumsummat(
    ::Type{TR}, n::Integer, x_min::TR, x_max::TR
) where {TR<:AbstractFloat}
    Q = cheb1_cumsummat(TR, n)
    Q .*= (x_max - x_min) / 2
    return Q
end

function cheb1_cumsummat(n::Integer, x_min::Float64, x_max::Float64)
    return cheb1_cumsummat(Float64, n, x_min, x_max)
end

export cheb1_cumsummat
