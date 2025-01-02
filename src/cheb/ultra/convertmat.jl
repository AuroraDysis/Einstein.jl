"""
    ultra_spconvert([TR=Float64], n::TI, λ::TR) where {TR<:AbstractFloat,TI<:Integer,TR<:Real}

Compute sparse representation for conversion operators.
Returns the truncation of the operator that transforms \$C^{\\lambda}\$
(Ultraspherical polynomials) to \$C^{\\lambda + 1}\$. The truncation gives
back a matrix of size n x n.

# References
- [chebfun/@ultraS/spconvert.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40ultraS/spconvert.m)
"""
function ultra_spconvert(::Type{TR}, n::TI, λ::TI) where {TI<:Integer,TR<:AbstractFloat}
    λf = convert(TR, λ)
    if λ == 0
        half = one(TR) / 2
        dg = fill(half, n - 2)
        # Construct sparse matrix with diagonals at 0 and 2
        return spdiagm(0 => [one(TR); half; dg], 2 => -dg)
    else
        dg = [λf / (λf + i) for i in 2:(n - 1)]
        return spdiagm(0 => [one(TR); λf / (λf + 1); dg], 2 => -dg)
    end
end

function ultra_spconvert(n::TI, λ::TR) where {TR<:Real,TI<:Integer}
    return ultra_spconvert(Float64, n, λ)
end

"""
    ultra_convertmat([TR=Float64], n::TI, K1::TI, K2::TI) where {TR<:AbstractFloat,TI<:Integer}

Conversion matrix used in the ultraspherical spectral method.
Returns N-by-N matrix realization of conversion operator between ultraspherical polynomial bases.
Maps N coefficients from C^{(K1)} basis to C^{(K2)} basis.

# References
- [chebfun/@ultraS/convertmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40ultraS/convertmat.m)
"""
function ultra_convertmat(
    ::Type{TR}, n::TI, K1::TI, K2::TI
) where {TR<:AbstractFloat,TI<:Integer}
    @argcheck n >= 2 "n must be positive"
    @argcheck K2 >= K1 "K2 must be greater than or equal to K1"
    S = sparse(I, n, n)
    for s in K1:(K2 - 1)
        S = ultra_spconvert(TR, n, s) * S
    end
    return S
end

function ultra_convertmat(n::TI, K1::TR, K2::TR) where {TI<:Integer,TR<:Real}
    return ultra_convertmat(Float64, n, K1, K2)
end

export ultra_convertmat, ultra_spconvert

@testset "ultra_convertmat" begin
    for type in [Float64, BigFloat]
        tol = 10 * eps(type)
        for n in 2:10
            for K1 in 1:3
                for K2 in K1:3
                    dom = -one(type) .. one(type)
                    @test isapprox(
                        ultra_convertmat(type, n, K1, K2),
                        Conversion(Ultraspherical(K1, dom), Ultraspherical(K2, dom))[
                            1:n, 1:n
                        ],
                        atol=tol,
                    )
                end
            end
        end
    end
end
