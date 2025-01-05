"""
MIT License

Copyright (c) 2025 Zhen Zhong
Copyright (c) 2019 Leo C. Stein

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

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

@doc raw"""
    sws_A0(::Type{TR}, s::Integer, l::Integer) where {TR<:AbstractFloat}

Calculate angular separation constant at a=0.
```math
{}_s A_{\ell m}(0)=l(l+1)-s(s+1)
```

# Arguments
- `TR`: type for floating point conversion
- `s::Integer`: spin
- `l::Integer`: angular number
"""
@inline function sws_A0(::Type{TR}, s::Integer, l::Integer) where {TR<:AbstractFloat}
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
        return sws_A0(TR, s, lprime) - c^2 * sws_calB(TR, s, lprime, m) +
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
    sws_eigM!(M::AbstractMatrix{TR}, s::Integer, c::Complex, m::Integer, l_max::Integer) where {TR<:AbstractFloat}

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
    return sws_eigM!(M, s, c, m, l_max)
end

function sws_eigM!(
    M::AbstractMatrix{Complex{TR}}, s::Integer, c::Complex{TR}, m::Integer, l_max::Integer
) where {TR<:AbstractFloat}
    l_range = sws_l_range(s, m, l_max)
    @argcheck size(M) == (length(l_range), length(l_range))
    @inbounds for (i, l) in enumerate(l_range), (j, lprime) in enumerate(l_range)
        M[i, j] = sws_Melem(TR, s, c, m, l, lprime)
    end
    return M
end

function sws_eigvals(
    ::Type{TR}, s::Integer, c::Complex{TR}, m::Integer, l_max::Integer
) where {TR<:AbstractFloat}
    M = sws_eigM(TR, s, c, m, l_max)
    evals = eigvals!(M; sortby=abs)
    return evals
end

@inline function sws_eigvalidx(s::Integer, l::Integer, m::Integer)
    l_min = sws_l_min(s, m)
    return l - l_min + 1
end

@doc raw"""
    sws_eigen(::Type{TR}, s::Integer, c::Complex{TR}, m::Integer, l_max::Integer) where {TR<:AbstractFloat}

Calculate eigenvalues and eigenvectors of the spherical-spheroidal decomposition matrix.

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

@doc raw"""
    SWSFun{TR<:AbstractFloat}()
    (swsf::SWSFun{TR})(θ::TR)
    coefficients(swsf::SWSFun)

Spherical-weighted spheroidal harmonics.
The eigenvectors contain the $C$ coefficients in the equation:
```math
{}_s S_{\ell m}(x; c)=\sum_{\ell^{\prime}=\ell_{\min }}^{\infty} C_{\ell^{\prime} \ell m}(c) {}_s S_{\ell^{\prime} m}(x; 0)
```
where $C$ is normalized by
```math
\sum_{\ell^{\prime}=\ell_{\text {min }}}^{\ell_{\max }}\left|C_{\ell^{\prime} \ell m}(c)\right|^2=1
```
and the phase is chosen such that $C_{\ell^{\prime} \ell m}(c)$ is real for $\ell^{\prime}=\ell$ [Cook:2014cta](@cite).
The spin-weighted spheroidal harmonics are normalized such that
```math
\int_0^\pi {}_s S_{\ell m}(x; c) {}_s S^*_{\ell m}(x; c) \sin(\theta) d\theta = 1
```

# Arguments
- `s::Integer`: spin
- `c::Complex{TR}`: oblateness parameter
- `m::Integer`: azimuthal number
- `l::Integer`: angular number
- `l_max::Integer`: maximum angular number
"""
struct SWSFun{TR<:AbstractFloat}
    s::Integer
    c::Complex{TR}
    m::Integer
    l_min::Integer
    l_max::Integer
    coeffs::Vector{Complex{TR}}
    Y_idx::Vector{Int}
    sqrt_2pi::TR
    Y_storage

    function SWSFun{TR}(
        s::Integer, c::Complex{TR}, m::Integer, l::Integer, l_max::Integer
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
        coeffs = vecs[:, vals_idx]
        scale = abs(coeffs[l - l_min + 1]) / coeffs[l - l_min + 1]
        coeffs .*= scale
        sqrt_2pi = sqrt(2 * convert(TR, π))
        return new{TR}(s, c, m, l_min, l_max, coeffs, Y_idx, sqrt_2pi, Y_storage)
    end
end

function (f::SWSFun{TR})(θ::TR) where {TR<:AbstractFloat}
    Y = sYlm_values!(f.Y_storage, θ, zero(TR), f.s)
    return f.sqrt_2pi * dot(f.coeffs, @view(Y[f.Y_idx]))
end

function coefficients(f::SWSFun)
    return f.coeffs
end

export sws_l_min, sws_A0, sws_eigM, sws_eigM!, sws_eigvals, sws_eigvalidx, sws_eigen, SWSFun
