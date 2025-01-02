"""
    spconvert([TR=Float64], n::TI, λ::TR) where {TR<:AbstractFloat,TI<:Integer,TR<:Real}

Compute sparse representation for conversion operators.
Returns the truncation of the operator that transforms \$C^{\\lambda}\$
(Ultraspherical polynomials) to \$C^{\\lambda + 1}\$. The truncation gives
back a matrix of size n x n.
"""
function spconvert(::Type{TR}, n::TI, λ::TI) where {TI<:Integer,TR<:AbstractFloat}
    λf = convert(TR, λ)
    if λ == 0
        half = one(TR) / 2
        dg = fill(half, n - 2)
        # Construct sparse matrix with diagonals at 0 and 2
        return spdiagm(0 => [one(TR); half; dg], 2 => [zero(TR); zero(TR); -dg])
    else
        dg = λf ./ (λf .+ (2:(n - 1)))
        return spdiagm(0 => [one(TR); λf / (λf + 1); dg], 2 => [zero(TR); zero(TR); -dg])
    end
end

function spconvert(n::TI, λ::TR) where {TR<:Real,TI<:Integer}
    return spconvert(Float64, n, λ)
end

"""
    convertmat([TR=Float64], n::TI, K1::TI, K2::TI) where {TR<:AbstractFloat,TI<:Integer}

Conversion matrix used in the ultraspherical spectral method.
Returns N-by-N matrix realization of conversion operator between ultraspherical polynomial bases.
Maps N coefficients from C^{(K1)} basis to C^{(K2 + 1)} basis.
Returns identity matrix if K2 < K1.
"""
function convertmat(::Type{TR}, n::TI, K1::TI, K2::TI) where {TR<:AbstractFloat,TI<:Integer}
    S = sparse(I, n, n)
    if K2 >= K1
        for s in K1:K2
            S = spconvert(TR, n, s) * S
        end
    end
    return S
end

function convertmat(n::TI, K1::TR, K2::TR) where {TI<:Integer,TR<:Real}
    return convertmat(Float64, n, K1, K2)
end

export convertmat, spconvert
