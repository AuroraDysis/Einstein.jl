using SphericalFunctions

"""
    sws_l_min(s::Integer, m::Integer)

Minimum allowed value of l for given s, m. Returns max(|s|, |m|).

# Arguments
- `s::Integer`: spin
- `m::Integer`: azimuthal number
"""
@inline function sws_l_min(s::Integer, m::Integer)
    return max(abs(s), abs(m))
end

"""
    sws_A0(::Type{TR}, s::Integer, l::Integer, m::Integer) where {TR<:AbstractFloat}

Calculate angular separation constant at a=0. Formula is A₀ = l(l+1) - s(s+1).

# Arguments
- `TR`: type for floating point conversion
- `s::Integer`: spin
- `l::Integer`: angular number
- `m::Integer`: azimuthal number
"""
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

@inline sws_l_range(s::Integer, m::Integer, l_max::Integer) = sws_l_min(s, m):l_max

"""
    sws_eigM(::Type{TR}, s::Integer, c::Complex{TR}, m::Integer, l_max::Integer) where {TR<:AbstractFloat}

Construct the spherical-spheroidal decomposition matrix truncated at l_max.

# Arguments
- `TR`: Type for floating point conversion
- `s::Integer`: spin
- `c::Complex`: oblateness parameter
- `m::Integer`: azimuthal number
- `l_max::Integer`: maximum angular number

# References
- [Cook:2014cta](@citet*)
- [duetosymmetry/qnm](https://github.com/duetosymmetry/qnm)
"""
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

@inline function sws_eigvalidx(s::Integer, l::Integer, m::Integer)
    l_min = sws_l_min(s, m)
    return l - l_min + 1
end

@doc raw"""
    sws_eigen(::Type{TR}, s::Integer, c::Complex{TR}, m::Integer, l_max::Integer) where {TR<:AbstractFloat}

Calculate eigenvalues and eigenvectors of the spherical-spheroidal decomposition matrix.
The eigenvectors contain the $C$ coefficients in the equation:
```math
{}_s S_{\ell m}(x; c)=\sum_{\ell^{\prime}=\ell_{\min }}^{\infty} C_{\ell^{\prime} \ell m}(c) {}_s S_{\ell^{\prime} m}(x; 0)
```
where $C$ is normalized by
```math
\sum_{\ell^{\prime}=\ell_{\text {min }}}^{\ell_{\max }}\left|C_{\ell^{\prime} \ell m}(c)\right|^2=1
```
and the phase is chosen such that $C_{\ell^{\prime} \ell m}(c)$ is real for $\ell^{\prime}=\ell$ [Cook:2014cta](@cite).

# Arguments
- `TR`: Type for floating point conversion
- `s::Integer`: spin
- `c::Complex`: oblateness parameter
- `m::Integer`: azimuthal number
- `l_max::Integer`: maximum angular number

# References
- [Cook:2014cta](@citet*)
- [duetosymmetry/qnm](https://github.com/duetosymmetry/qnm)
"""
function sws_eigen(
    ::Type{TR}, s::Integer, c::Complex{TR}, m::Integer, l_max::Integer
) where {TR<:AbstractFloat}
    M = sws_eigM(TR, s, c, m, l_max)
    vals, vecs = eigen!(M)
    return vals, vecs
end

struct SWSFun{TR<:AbstractFloat}
    s::Integer
    c::Complex{TR}
    m::Integer
    l_min::Integer
    l_max::Integer
    coeffs::Vector{Complex{TR}}
    Y_idx::Vector{Int}
    Y_storage

    function SWSFun(
        ::Type{TR}, s::Integer, c::Complex{TR}, m::Integer, l::Integer, l_max::Integer
    ) where {TR<:AbstractFloat}
        l_min = sws_l_min(s, m)
        Y_storage = sYlm_prep(l_max, s, TR, l_min)
        vals, vecs = sws_eigen(TR, s, c, m, l_max)
        vals_idx = sws_eigvalidx(s, l, m)
        i = 1
        Y_idx = zeros(Int, l_max - l_min + 1)
        # TODO: make this more efficient
        for lprime in l_min:l_max
            for mprime in (-lprime):lprime
                if m == mprime
                    Y_idx[lprime - l_min + 1] = i
                end
                i += 1
            end
        end
        return new{TR}(s, c, m, l_min, l_max, vecs[:, vals_idx], Y_idx, Y_storage)
    end
end

function (f::SWSFun)(θ::TR) where {TR<:AbstractFloat}
    Y = sYlm_values!(f.Y_storage, θ, zero(TR), f.s)
    return dot(f.coeffs, @view(Y[f.Y_idx]))
end

export sws_l_min, sws_A0, sws_eigM, sws_eigvals, sws_eigvalidx, sws_eigen, SWSFun
