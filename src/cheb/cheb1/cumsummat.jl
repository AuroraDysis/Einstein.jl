function cheb1_cumsummat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    A = cheb1_amat(TR, n)
    S = cheb1_smat(TR, n)
    B = cheb_coeffs_cumsummat(TR, n)
    return S * B * A
end

function cheb1_cumsummat(
    ::Type{TR}, n::TI, x_min::TR, x_max::TR
) where {TR<:AbstractFloat,TI<:Integer}
    A = cheb1_amat(TR, n)
    S = cheb1_smat(TR, n)
    B = cheb_coeffs_cumsummat(TR, n, x_min, x_max)
    return S * B * A
end

export cheb1_cumsummat
