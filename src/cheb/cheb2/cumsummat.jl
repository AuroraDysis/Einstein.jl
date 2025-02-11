"""
    cheb2_cumsummat([TF=Float64], n::Integer) where {TF<:AbstractFloat}
    cheb2_cumsummat([TF=Float64], n::Integer, x_min::TF, x_max::TF) where {TF<:AbstractFloat}

Compute Chebyshev integration matrix that maps function values
at `n` Chebyshev points of the 2st kind to values of the integral of the interpolating
polynomial at those points, with the convention that the first value is zero.

# References
- [chebfun/@chebcolloc2/chebcolloc2.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc2/chebcolloc2.m)
"""
function cheb2_cumsummat(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    A = cheb2_analysis_matrix(TF, n)
    S = cheb2_synthesis_matrix(TF, n)
    B = cheb_coeffs_cumsummat(TF, n)
    Q = S * B * A
    @inbounds Q[1, :] .= 0
    return Q
end

function cheb2_cumsummat(n::Integer)
    return cheb2_cumsummat(Float64, n)
end

function cheb2_cumsummat(::Type{TF}, n::Integer, x_min::TF, x_max::TF) where {TF<:AbstractFloat}
    Q = cheb2_cumsummat(TF, n)
    scale = (x_max - x_min) / 2
    Q .*= scale
    return Q
end

function cheb2_cumsummat(n::Integer, x_min::Float64, x_max::Float64)
    return cheb2_cumsummat(Float64, n, x_min, x_max)
end

export cheb2_cumsummat
