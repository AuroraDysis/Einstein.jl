function cheb2_cumsummat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    A = cheb2_amat(TR, n)
    S = cheb2_smat(TR, n)
    B = cheb_coeffs_cumsummat(TR, n)
    @inbounds B[1, :] .= 0
    return S * B * A
end

function cheb2_cumsummat(
    ::Type{TR}, n::TI, x_min::TR, x_max::TR
) where {TR<:AbstractFloat,TI<:Integer}
    A = cheb2_amat(TR, n)
    S = cheb2_smat(TR, n)
    B = cheb_coeffs_cumsummat(TR, n, x_min, x_max)
    @inbounds B[1, :] .= 0
    return S * B * A
end

export cheb2_cumsummat
