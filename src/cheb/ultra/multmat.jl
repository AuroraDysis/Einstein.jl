"""
    ultra_multmat(a::VT, λ::TI) where {TR<:AbstractFloat,VT<:AbstractVector{TR},TI<:Integer}

Construct nxn multiplication matrix
representing the multiplication of F in the \$C^{(\\lambda)}\$ basis.

# Arguments
- `a::VT` : Vector of Chebyshev coefficients (!a is modified in place!)
- `λ::TI` : Order of the ultraspherical basis

# References
- [chebfun/@ultraS/multmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40ultraS/multmat.m)
"""
function ultra_multmat(
    a::VT, λ::TI
) where {TR<:AbstractFloat,VT<:AbstractVector{TR},TI<:Integer}
    n = length(a)

    if λ == 0 # λ == 0 case (Chebyshev T)
        half = one(TR) / 2
        @inbounds a[2:end] .*= half
        M = sparse(Toeplitz(a, a))
        H = ultra_sphankel(@view(a[2:end]))
        M[2:n, 1:(n - 1)] .+= H
        return M
    elseif λ == 1 # λ == 1 case (Chebyshev U)
        half = one(TR) / 2
        @inbounds a[2:end] .*= half
        M = sparse(Toeplitz(a, a) .* half)
        M[1:(n - 2), 1:(n - 2)] .-= ultra_sphankel(@view(a[3:end]))
        return M
    else # General λ case, Convert from T to C^(λ-1)
        a = ultra_convertmat(TR, n, 0, λ) * a
        m = 2n
        M0 = sparse(I, m, m)

        d1 = ((2λ):(2λ + m - 2)) ./ (2((λ + 1):(λ + m - 1)))
        d2 = (1:(m - 1)) ./ (2(λ:(λ + m - 2)))

        Mx = spdiagm(-1 => d2, 1 => d1)
        M1 = 2λ .* Mx

        M = a[1] * M0 + a[2] * M1
        for nn in 1:(n - 2)
            M2 = 2(nn + λ) / (nn + 1) * Mx * M1 .- (nn + 2λ - 1) / (nn + 1) * M0
            M .+= a[nn + 2] * M2
            M0, M1 = M1, M2
        end
        return M[1:n, 1:n]
    end
end
