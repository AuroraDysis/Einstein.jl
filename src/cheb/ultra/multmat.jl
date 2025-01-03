"""
    ultra_multmat(coeffs::VT, λ::TI) where {TR<:AbstractFloat,VT<:AbstractVector{TR},TI<:Integer}

Construct nxn multiplication matrix
representing the multiplication of F in the \$C^{(\\lambda)}\$ basis.

# Arguments
- `coeffs::VT` : Vector of Chebyshev coefficients
- `λ::TI` : Order of the ultraspherical basis

# References
- [chebfun/@ultraS/multmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40ultraS/multmat.m)
"""
function ultra_multmat(
    coeffs::VT, λ::TI
) where {TR<:AbstractFloat,VT<:AbstractVector{TR},TI<:Integer}
    n = length(coeffs)

    coeffs = deepcopy(coeffs)

    if λ == 0 # λ == 0 case (Chebyshev T)
        half = one(TR) / 2
        @inbounds coeffs[2:end] .*= half
        M = Matrix(Toeplitz(coeffs, coeffs))
        H = ultra_sphankel(@view(coeffs[2:end]))
        M[2:n, 1:(n - 1)] .+= H
        return M
    elseif λ == 1 # λ == 1 case (Chebyshev U)
        half = one(TR) / 2
        @inbounds coeffs[2:end] .*= half
        M = Matrix(Toeplitz(coeffs, coeffs))
        M[1:(n - 2), 1:(n - 2)] .-= ultra_sphankel(@view(coeffs[3:end]))
        return M
    else # General λ case, Convert from T to C^λ
        coeffs = ultra_convertmat(TR, n, 0, λ) * coeffs
        m = 2n
        M0 = BandedMatrix(Eye{TR}(m))

        λf = convert(TR, λ)
        d1 = Array{TR}(undef, m - 1)
        d2 = Array{TR}(undef, m - 1)
        @inbounds for i in 0:(m - 2)
            d1[i + 1] = (2λf + i) / (2(λf + i + 1))
            d2[i + 1] = (1 + i) / (2(λf + i))
        end

        Mx = BandedMatrix(-1 => d2, 1 => d1)
        M1 = 2λ .* Mx

        M = BandedMatrix{TR}(undef, (m, m), (n, n))
        @. M = coeffs[1] * M0 + coeffs[2] * M1
        for nn in 1:(n - 2)
            M2 = 2(nn + λf) / (nn + 1) * Mx * M1 .- (nn + 2λf - 1) / (nn + 1) * M0
            @. M += coeffs[nn + 2] * M2
            M0, M1 = M1, M2
        end

        return Matrix(@view(M[1:n, 1:n]))
    end
end

export ultra_multmat
