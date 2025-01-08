@doc raw"""
    ultra_multmat(coeffs::AbstractVector{TRC}, λ::Integer) where {TRC<:Union{AbstractFloat,Complex{<:AbstractFloat}}
    ultra_multmat!(coeffs::AbstractVector{TRC}, λ::Integer) where {TRC<:Union{AbstractFloat,Complex{<:AbstractFloat}}

Construct nxn multiplication matrix representing the multiplication of F in the $C^{(\lambda)}$ basis.
`ultra_multmat!` will overwrite the input vector `coeffs`.

# Arguments
- `coeffs::AbstractVector{TRC}` : Vector of Chebyshev coefficients
- `λ::Integer` : Order of the ultraspherical basis

# References
- [chebfun/@ultraS/multmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40ultraS/multmat.m)
"""
function ultra_multmat!(
    coeffs::AbstractVector{TRC}, λ::Integer
) where {TRC<:AbstractFloatOrComplex}
    n = length(coeffs)

    if λ == 0 # λ == 0 case (Chebyshev T)
        half = one(real(TRC)) / 2
        @inbounds coeffs[2:end] .*= half
        M = Matrix(Toeplitz(coeffs, coeffs))
        H = ultra_sphankel(@view(coeffs[2:end]))
        M[2:n, 1:(n - 1)] .+= H
        return M
    elseif λ == 1 # λ == 1 case (Chebyshev U)
        half = one(real(TRC)) / 2
        @inbounds coeffs[2:end] .*= half
        M = Matrix(Toeplitz(coeffs, coeffs))
        M[1:(n - 2), 1:(n - 2)] .-= ultra_sphankel(@view(coeffs[3:end]))
        return M
    else # General λ case, Convert from T to C^λ
        coeffs = ultra_convertmat(real(TRC), 0, λ, n) * coeffs
        m = 2n
        M0 = BandedMatrix(Eye{TRC}(m))

        λf = convert(TRC, λ)
        d1 = Array{TRC}(undef, m - 1)
        d2 = Array{TRC}(undef, m - 1)
        @inbounds for i in 0:(m - 2)
            d1[i + 1] = (2λf + i) / (2(λf + i + 1))
            d2[i + 1] = (1 + i) / (2(λf + i))
        end

        Mx = BandedMatrix(-1 => d2, 1 => d1)
        M1 = 2λ .* Mx

        M = BandedMatrix{TRC}(undef, (m, m), (n, n))
        @. M = coeffs[1] * M0 + coeffs[2] * M1
        for nn in 1:(n - 2)
            M2 = 2(nn + λf) / (nn + 1) * Mx * M1 .- (nn + 2λf - 1) / (nn + 1) * M0
            @. M += coeffs[nn + 2] * M2
            M0, M1 = M1, M2
        end

        return Matrix(@view(M[1:n, 1:n]))
    end
end

function ultra_multmat(
    coeffs::AbstractVector{TRC}, λ::Integer
) where {TRC<:AbstractFloatOrComplex}
    return ultra_multmat!(deepcopy(coeffs), λ)
end

export ultra_multmat, ultra_multmat!
