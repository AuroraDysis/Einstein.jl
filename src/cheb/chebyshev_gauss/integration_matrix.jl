"""
    cheb_gauss_integration_matrix([TF=Float64], n::Integer) where {TF<:AbstractFloat}
    cheb_gauss_integration_matrix([TF=Float64], n::Integer, domain_width::TF) where {TF<:AbstractFloat}

Compute Chebyshev integration matrix that maps function values
at `n` Chebyshev points of the 1st kind to values of the integral of the interpolating
polynomial at those points, with the convention that the first value is zero.

# References
- [chebfun/@chebcolloc1/chebcolloc1.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebcolloc1/chebcolloc1.m)
"""
function cheb_gauss_integration_matrix(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    A = cheb_gauss_vals2coeffs_matrix(TF, n)
    S = cheb_gauss_coeffs2vals_matrix(TF, n)
    B = cheb_series_integration_matrix(TF, n)
    Q = S * B * A
    return Q
end

function cheb_gauss_integration_matrix(n::Integer)
    return cheb_gauss_integration_matrix(Float64, n)
end

function cheb_gauss_integration_matrix(
    ::Type{TF}, n::Integer, domain_width::TF
) where {TF<:AbstractFloat}
    Q = cheb_gauss_integration_matrix(TF, n)
    jacobian = domain_width / 2
    Q .*= jacobian
    return Q
end

function cheb_gauss_integration_matrix(n::Integer, domain_width::Float64)
    return cheb_gauss_integration_matrix(Float64, n, domain_width)
end

export cheb_gauss_integration_matrix
