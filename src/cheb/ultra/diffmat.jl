@doc raw"""
    ultra_diffmat([TR=Float64], m::Integer, n::Integer) where TR<:AbstractFloat

Differentiation matrices for ultraspherical spectral method that takes $n$ Chebyshev
coefficients and returns $n$ $C^{(m)}$ coefficients that represent the derivative
of the Chebyshev series. Here, $C^{(k)}$ is the ultraspherical polynomial basis
with parameter $k$.

# Arguments
- `m::Integer`: Order of differentiation
- `n::Integer`: Number of points

# References
- [chebfun/@ultraS/diffmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40ultraS/diffmat.m)
"""
function ultra_diffmat(::Type{TR}, m::TI, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    @argcheck m >= 1 "m must be nonnegative"
    @argcheck n >= 2 "n must be positive"

    nm1 = n - 1
    D = spdiagm(1 => collect(TR, 1:nm1))
    for s in 1:(m - 1)
        D = spdiagm(1 => 2s * ones(TR, nm1)) * D
    end

    return D
end

function ultra_diffmat(m::TI, n::TI) where {TI<:Integer}
    return ultra_diffmat(Float64, m, n)
end

function ultra_diffmat(
    ::Type{TR}, m::TI, n::TI, lower_bound::TR, upper_bound::TR
) where {TR<:AbstractFloat,TI<:Integer}
    @argcheck lower_bound < upper_bound "lower_bound must be less than upper_bound"
    @argcheck n >= 2 "n must be positive"

    scale = (2 / (upper_bound - lower_bound))^m
    D = ultra_diffmat(TR, m, n)
    D .*= scale

    return D
end

function ultra_diffmat(m::TI, n::TI, lower_bound::Float64, upper_bound::Float64) where {TI<:Integer}
    return ultra_diffmat(Float64, m, n, lower_bound, upper_bound)
end

export ultra_diffmat
