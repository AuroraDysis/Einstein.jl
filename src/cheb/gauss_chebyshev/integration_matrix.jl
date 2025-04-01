"""
    gauss_chebyshev_integration_matrix([TF=Float64], n::Integer) where {TF<:AbstractFloat}
    gauss_chebyshev_integration_matrix([TF=Float64], n::Integer, lower_bound::TF, upper_bound::TF) where {TF<:AbstractFloat}

Compute Chebyshev integration matrix that maps function values
at `n` Chebyshev points of the 1st kind to values of the integral of the interpolating
polynomial at those points, with the convention that the first value is zero.

# References
- [chebfun/@chebcolloc1/chebcolloc1.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc1/chebcolloc1.m)
"""
function gauss_chebyshev_integration_matrix(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    A = gauss_chebyshev_vals2coeffs_matrix(TF, n)
    S = gauss_chebyshev_coeffs2vals_matrix(TF, n)
    B = chebyshevt_integration_matrix(TF, n)
    Q = S * B * A
    return Q
end

function gauss_chebyshev_integration_matrix(n::Integer)
    return gauss_chebyshev_integration_matrix(Float64, n)
end

function gauss_chebyshev_integration_matrix(
    ::Type{TF}, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    Q = gauss_chebyshev_integration_matrix(TF, n)
    Q .*= (upper_bound - lower_bound) / 2
    return Q
end

function gauss_chebyshev_integration_matrix(n::Integer, lower_bound::Float64, upper_bound::Float64)
    return gauss_chebyshev_integration_matrix(Float64, n, lower_bound, upper_bound)
end

export gauss_chebyshev_integration_matrix
