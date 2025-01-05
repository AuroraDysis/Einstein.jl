@inline function sws_l_min(s::T, m::T) where {T<:Integer}
    return max(abs(s), abs(m))
end

@inline function sws_A0(
    ::Type{TR}, s::Integer, l::Integer, m::Integer
) where {TR<:AbstractFloat}
    return convert(TR, l * (l + 1) - s * (s + 1))
end

@inline function sws_recF(
    ::Type{TR}, s::Integer, l::Integer, m::Integer
) where {TR<:AbstractFloat}
    if (s == 0 && l + 1 == 0)
        return zero(TR)
    end

    num1 = convert(TR, (l + 1)^2 - m^2)
    den1 = convert(TR, (2l + 3) * (2l + 1))

    num2 = convert(TR, (l + 1)^2 - s^2)
    den2 = convert(TR, (l + 1)^2)

    return sqrt(num1 / den1) * sqrt(num2 / den2)
end

@inline function sws_recG(
    ::Type{TR}, s::Integer, l::Integer, m::Integer
) where {TR<:AbstractFloat}
    if l == 0
        return zero(TR)
    end

    num1 = convert(TR, l^2 - m^2)
    den1 = convert(TR, (2l - 1) * (2l + 1))

    num2 = convert(TR, l^2 - s^2)
    den2 = convert(TR, l^2)

    return sqrt(num1 / den1) * sqrt(num2 / den2)
end

@inline function sws_recH(
    ::Type{TR}, s::Integer, l::Integer, m::Integer
) where {TR<:AbstractFloat}
    if l == 0 || s == 0
        return zero(TR)
    end

    num = convert(TR, m * s)
    den = convert(TR, l * (l + 1))

    return -num / den
end

@inline function sws_calA(
    ::Type{TR}, s::Integer, l::Integer, m::Integer
) where {TR<:AbstractFloat}
    return sws_recF(TR, s, l, m) * sws_recF(TR, s, l + 1, m)
end

@inline function sws_calD(
    ::Type{TR}, s::Integer, l::Integer, m::Integer
) where {TR<:AbstractFloat}
    return sws_recF(TR, s, l, m) * (sws_recH(TR, s, l + 1, m) + sws_recH(TR, s, l, m))
end

@inline function sws_calB(
    ::Type{TR}, s::Integer, l::Integer, m::Integer
) where {TR<:AbstractFloat}
    return sws_recF(TR, s, l, m) * sws_recG(TR, s, l + 1, m) +
           sws_recG(TR, s, l, m) * sws_recF(TR, s, l - 1, m) +
           sws_recH(TR, s, l, m)^2
end

@inline function sws_calE(
    ::Type{TR}, s::Integer, l::Integer, m::Integer
) where {TR<:AbstractFloat}
    return sws_recG(TR, s, l, m) * (sws_recH(TR, s, l - 1, m) + sws_recH(TR, s, l, m))
end

@inline function sws_calC(
    ::Type{TR}, s::Integer, l::Integer, m::Integer
) where {TR<:AbstractFloat}
    return sws_recG(TR, s, l, m) * sws_recG(TR, s, l - 1, m)
end

function sws_Melem(
    ::Type{TR}, s::Integer, c::Complex{TR}, m::Integer, l::Integer, lprime::Integer
) where {TR<:AbstractFloat}
    if lprime == l - 2
        return -c^2 * sws_calA(TR, s, lprime, m)
    elseif lprime == l - 1
        return -c^2 * sws_calD(TR, s, lprime, m) + 2 * c * s * sws_recF(TR, s, lprime, m)
    elseif lprime == l
        return sws_A0(TR, s, lprime, m) - c^2 * sws_calB(TR, s, lprime, m) +
               2 * c * s * sws_recH(TR, s, lprime, m)
    elseif lprime == l + 1
        return -c^2 * sws_calE(TR, s, lprime, m) + 2 * c * s * sws_recG(TR, s, lprime, m)
    elseif lprime == l + 2
        return -c^2 * sws_calC(TR, s, lprime, m)
    end

    return zero(Complex{TR})
end

@inline sws_l_range(s::Integer, m::Integer, l_max::Integer) = sws_l_min(s, m):(l_max + 1)

function sws_eigM(
    ::Type{TR}, s::Integer, c::Complex{TR}, m::Integer, l_max::Integer
) where {TR<:AbstractFloat}
    l_range = sws_l_range(s, m, l_max)
    M_size = length(l_range)
    M = Array{Complex{TR}}(undef, M_size, M_size)
    for (i, l) in enumerate(l_range), (j, lprime) in enumerate(l_range)
        M[i, j] = sws_Melem(TR, s, c, m, l, lprime)
    end
    return M
end

function sws_eigvals(
    ::Type{TR}, s::Integer, c::Complex{TR}, m::Integer, l_max::Integer
) where {TR<:AbstractFloat}
    M = sws_eigM(TR, s, c, m, l_max)
    vals = eigvals!(M)
    return vals
end

@inline function sws_eigvalidx(s::Integer, m::Integer, l::Integer)
    l_min = sws_l_min(s, m)
    return l - l_min + 1
end

function sws_eigen(
    ::Type{TR}, s::Integer, c::Complex{TR}, m::Integer, l_max::Integer
) where {TR<:AbstractFloat}
    M = sws_eigM(TR, s, c, m, l_max)
    vals, vecs = eigen!(M)
    return vals, vecs
end

export sws_l_min, sws_A0, sws_eigM, sws_eigvals, sws_eigvalidx, sws_eigen
