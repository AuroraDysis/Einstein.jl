"""
    ultra_multmat(a::VT, λ::TI) where {TR<:AbstractFloat,VT<:AbstractVector{TR},TI<:Integer}

Construct nxn multiplication matrix
representing the multiplication of F in the \$C^{(\\lambda)}\$ basis.

# References
- [chebfun/@ultraS/multmat.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40ultraS/multmat.m)
"""
function ultra_multmat(
    a::VT, λ::TI
) where {TR<:AbstractFloat,VT<:AbstractVector{TR},TI<:Integer}
    n = length(a)

    if λ == 0 # λ == 0 case (Chebyshev T)
        a .*= 0.5
        M = sparse(Toeplitz([2a[1]; a[2:end]], [2a[1]; a[2:end]]))
        H = ultra_sphankel(a[2:end])
        sub1 = 2:n
        sub2 = 1:(n - 1)
        M[sub1, sub2] .+= H
        return M
    elseif λ == 1 # λ == 1 case (Chebyshev U)
        M = sparse(Toeplitz([2a[1]; a[2:end]], [2a[1]; a[2:end]]) ./ 2)
        sub = 1:(n - 2)
        M[sub, sub] .-= ultra_sphankel(a[3:end] ./ 2)
        return M
    else # General λ case, Convert from T to C^(λ-1)
        a = ultra_convertmat(n, 0, λ) * a
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
