"""
    cheb2_cumsummat([T=Float64], n::Integer) where {T<:AbstractFloat}
    cheb2_cumsummat([T=Float64], n::Integer, x_min::T, x_max::T) where {T<:AbstractFloat}

Compute Chebyshev integration matrix that maps function values
at `n` Chebyshev points of the 2st kind to values of the integral of the interpolating
polynomial at those points, with the convention that the first value is zero.

# References
- [chebfun/@chebcolloc2/chebcolloc2.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc2/chebcolloc2.m)
"""
function cheb2_cumsummat(::Type{T}, n::Integer) where {T<:AbstractFloat}
    A = cheb2_amat(T, n)
    S = cheb2_smat(T, n)
    B = cheb_coeffs_cumsummat(T, n)
    Q = S * B * A
    @inbounds Q[1, :] .= 0
    return Q
end

function cheb2_cumsummat(n::Integer)
    return cheb2_cumsummat(Float64, n)
end

function cheb2_cumsummat(n::Integer, x_min::T, x_max::T) where {T<:AbstractFloat}
    Q = cheb2_cumsummat(T, n)
    Q .*= (x_max - x_min) / 2
    return Q
end

function cheb2_cumsummat(::Type{T}, n::Integer, x_min::T, x_max::T) where {T<:AbstractFloat}
    return cheb2_cumsummat(n, x_min, x_max)
end

export cheb2_cumsummat
