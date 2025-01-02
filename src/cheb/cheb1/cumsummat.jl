function cheb1_cumsummat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    A = cheb1_amat(TR, n)
    S = cheb1_smat(TR, n)
    B = cheb_coeffs_cumsummat(TR, n)
    Q = S * B * A
    return Q
end

function cheb1_cumsummat(n::TI) where {TI<:Integer}
    return cheb1_cumsummat(Float64, n)
end

function cheb1_cumsummat(
    ::Type{TR}, n::TI, x_min::TR, x_max::TR
) where {TR<:AbstractFloat,TI<:Integer}
    Q = cheb1_cumsummat(TR, n)
    Q .*= (x_max - x_min) / 2
    return Q
end

function cheb1_cumsummat(n::TI, x_min::Float64, x_max::Float64) where {TI<:Integer}
    return cheb1_cumsummat(Float64, n, x_min, x_max)
end

export cheb1_cumsummat
