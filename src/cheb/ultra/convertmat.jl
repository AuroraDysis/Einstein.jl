@doc raw"""
    ultra_spconvert([T=Float64], n::Integer, λ::T) where {T<:AbstractFloat}

Compute sparse representation for conversion operators.
Returns the truncation of the operator that transforms $C^{\lambda}$
(Ultraspherical polynomials) to $C^{\lambda + 1}$. The truncation gives
back a matrix of size n x n.

# References
- [chebfun/@ultraS/spconvert.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40ultraS/spconvert.m)
"""
function ultra_spconvert(::Type{T}, n::Integer, λ::Integer) where {T<:AbstractFloat}
    λf = convert(T, λ)
    if λ == 0
        half = one(T) / 2
        dg = fill(half, n - 2)
        # Construct sparse matrix with diagonals at 0 and 2
        return spdiagm(0 => [one(T); half; dg], 2 => -dg)
    else
        dg = [λf / (λf + i) for i in 2:(n - 1)]
        return spdiagm(0 => [one(T); λf / (λf + 1); dg], 2 => -dg)
    end
end

function ultra_spconvert(n::Integer, λ::T) where {T<:Real,Integer<:Integer}
    return ultra_spconvert(Float64, n, λ)
end

"""
    ultra_convertmat([T=Float64], n::Integer, K1::Integer, K2::Integer) where {T<:AbstractFloat}

Conversion matrix used in the ultraspherical spectral method.
Returns N-by-N matrix realization of conversion operator between ultraspherical polynomial bases.
Maps N coefficients from C^{(K1)} basis to C^{(K2)} basis.

# References
- [chebfun/@ultraS/convertmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40ultraS/convertmat.m)
"""
function ultra_convertmat(
    ::Type{T}, n::Integer, K1::Integer, K2::Integer
) where {T<:AbstractFloat}
    @argcheck n >= 2 "n must be positive"
    @argcheck K2 >= K1 "K2 must be greater than or equal to K1"
    S = sparse(I, n, n)
    for s in K1:(K2 - 1)
        S = ultra_spconvert(T, n, s) * S
    end
    return S
end

function ultra_convertmat(n::Integer, K1::T, K2::T) where {T<:Real}
    return ultra_convertmat(Float64, n, K1, K2)
end

export ultra_convertmat, ultra_spconvert
