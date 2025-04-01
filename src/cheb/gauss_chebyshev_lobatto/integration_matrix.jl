"""
    gauss_chebyshev_lobatto_integration_matrix([TF=Float64], n::Integer) where {TF<:AbstractFloat}
    gauss_chebyshev_lobatto_integration_matrix([TF=Float64], n::Integer, lower_bound::TF, upper_bound::TF) where {TF<:AbstractFloat}

Compute Chebyshev integration matrix that maps function values
at `n` Chebyshev points of the 2st kind to values of the integral of the interpolating
polynomial at those points, with the convention that the first value is zero.

# References
- [chebfun/@chebcolloc2/chebcolloc2.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc2/chebcolloc2.m)
"""
function gauss_chebyshev_lobatto_integration_matrix(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    A = gauss_chebyshev_lobatto_vals2coeffs_matrix(TF, n)
    S = gauss_chebyshev_lobatto_coeffs2vals_matrix(TF, n)
    B = chebyshevt_integration_matrix(TF, n)
    Q = S * B * A
    @inbounds Q[1, :] .= 0
    return Q
end

function gauss_chebyshev_lobatto_integration_matrix(n::Integer)
    return gauss_chebyshev_lobatto_integration_matrix(Float64, n)
end

function gauss_chebyshev_lobatto_integration_matrix(
    ::Type{TF}, n::Integer, lower_bound::TF, upper_bound::TF
) where {TF<:AbstractFloat}
    Q = gauss_chebyshev_lobatto_integration_matrix(TF, n)
    scale = (upper_bound - lower_bound) / 2
    Q .*= scale
    return Q
end

function gauss_chebyshev_lobatto_integration_matrix(n::Integer, lower_bound::Float64, upper_bound::Float64)
    return gauss_chebyshev_lobatto_integration_matrix(Float64, n, lower_bound, upper_bound)
end

export gauss_chebyshev_lobatto_integration_matrix
