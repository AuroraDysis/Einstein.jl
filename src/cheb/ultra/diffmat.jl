@doc raw"""
    ultra_diffmat(n::Integer, m::Integer)

Differentiation matrices for ultraspherical spectral method that takes $n$ Chebyshev
coefficients and returns $n$ $C^{(m)}$ coefficients that represent the derivative
of the Chebyshev series. Here, $C^{(k)}$ is the ultraspherical polynomial basis
with parameter $k$.

# Arguments
- `n::Integer`: Number of points
- `m::Integer`: Order of differentiation

# References
- [chebfun/@ultraS/diffmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40ultraS/diffmat.m)
"""
function ultra_diffmat(n::TI, m::TI) where {TI<:Integer}
    @argcheck n >= 2 "n must be positive"
    @argcheck m >= 1 "m must be nonnegative"

    nm1 = n - 1
    D = spdiagm(1 => 1:nm1)
    for s in 1:(m - 1)
        D = spdiagm(1 => 2s * ones(TI, nm1)) * D
    end

    return D
end

export ultra_diffmat
