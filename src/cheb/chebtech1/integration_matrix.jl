"""
    cheb1_integration_matrix([TF=Float64], n::Integer) where {TF<:AbstractFloat}
    cheb1_integration_matrix([TF=Float64], n::Integer, lower_bound::TF, upper_bound::TF) where {TF<:AbstractFloat}

Compute Chebyshev integration matrix that maps function values
at `n` Chebyshev points of the 1st kind to values of the integral of the interpolating
polynomial at those points, with the convention that the first value is zero.

# References
- [chebfun/@chebcolloc1/chebcolloc1.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc1/chebcolloc1.m)
"""
function cheb1_integration_matrix(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    A = chebtech1_analysis_matrix(TF, n)
    S = cheb1_synthesis_matrix(TF, n)
    B = cheb_coeffs_integration_matrix(TF, n)
    Q = S * B * A
    return Q
end

function cheb1_integration_matrix(n::Integer)
    return cheb1_integration_matrix(Float64, n)
end

function cheb1_integration_matrix(
    ::Type{TF}, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    Q = cheb1_integration_matrix(TF, n)
    Q .*= (upper_bound - lower_bound) / 2
    return Q
end

function cheb1_integration_matrix(n::Integer, lower_bound::Float64, upper_bound::Float64)
    return cheb1_integration_matrix(Float64, n, lower_bound, upper_bound)
end

function _cheb_integration_matrix(
    ::ChebyshevT, ::Type{TR}, n::Integer, lower_bound::TR, upper_bound::TR
) where {TR<:AbstractFloat}
    return cheb1_integration_matrix(TR, n, lower_bound, upper_bound)
end

export cheb1_integration_matrix
