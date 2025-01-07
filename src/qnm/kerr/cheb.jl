@with_kw struct QNMKerrRadialChebCache{TR<:AbstractFloat}
    a::TR # black hole spin
    s::Integer # spin weight of the field
    m::Integer # azimuthal mode number
    A0::AbstractMatrix{Complex{TR}}
    A1::AbstractMatrix{Complex{TR}}
    A2::AbstractMatrix{TR}
    An::AbstractMatrix{Complex{TR}}
    Lm::AbstractMatrix{Complex{TR}}
end

"""
    qnm_kerr_radial_cheb_cache(::Type{TR}, a, s, m, n; ρ_min=0, ρ_max=1, lo_bc=BCType.Natural, hi_bc=BCType.Natural)

Initialize the cache for Kerr QNM nonlinear eigenvalue problem using Ultraspherical spectral method.

# Arguments
- `TR`: Type parameter for floating-point precision
- `a`: Black hole spin parameter
- `s`: Spin weight of the field
- `m`: Azimuthal mode number
- `n`: Number of Chebyshev points
- `ρ_min`: Minimum compactified radial coordinate (default: 0)
- `ρ_max`: Maximum compactified radial coordinate (default: 1)
- `lo_bc`: Boundary condition at ρ_min (default: Natural)
- `hi_bc`: Boundary condition at ρ_max (default: Natural) - not implemented yet

# Returns
- `QNMKerrRadialChebCache` object containing pre-computed matrices

# References
- [Ripley:2022ypi](@citet*)
"""
function qnm_kerr_radial_cheb_cache(
    ::Type{TR},
    a::TR,
    s::Integer,
    m::Integer,
    n::Integer;
    ρ_min::TR=zero(TR),
    ρ_max::TR=one(TR),
    lo_bc::BCType.T=BCType.Natural,
    hi_bc::BCType.T=BCType.Natural,
) where {TR<:AbstractFloat}
    # use M = 1 unit
    M = one(TR)

    dom = ρ_min .. ρ_max
    chebSpace = Chebyshev(dom)
    ultraSpace = Ultraspherical(2, dom)
    conversion = Conversion(chebSpace, ultraSpace)
    conversionA1 = Conversion(Ultraspherical(1, dom), ultraSpace)

    ρ = Fun(chebSpace)
    c02 = -ρ^2 * (1 - 2 * M * ρ + a^2 * ρ^2)
    c01 = -2 * ρ * (1 + s - (1im * a * m + M * (3 + s)) * ρ + 2 * a^2 * ρ^2)
    c00 = 2 * ρ * (1im * a * m + M * (1 + s) - a^2 * ρ)
    A0c = (c02 * 𝒟^2 + c01 * 𝒟 + c00):chebSpace

    c11 = 2im * (1 + ρ^2 * (-8 * M^2 + a^2 * (1 + 4 * M * ρ)))
    c10 =
        2im * a^2 * ρ * (1 + 6 * M * ρ) +
        2 * a * m * (1 + 4 * M * ρ) +
        4im * M * (s - 2 * M * (s + 2) * ρ)
    A1 = (c11 * 𝒟 + c10):chebSpace
    A1c = conversionA1 * A1

    c20 = -16 * M^2 * (1 + 2 * M * ρ) + a^2 * (1 + 4 * M * ρ)^2
    A2 = (c20 * conversion):chebSpace

    A0m = Matrix(@view(A0c[1:n, 1:n]))
    A1m = Matrix(@view(A1c[1:n, 1:n]))
    A2m = Matrix(@view(A2[1:n, 1:n]))
    Anm = Matrix(@view(conversion[1:n, 1:n]))

    if lo_bc == BCType.Dirichlet
        A0m[end, :] .= TR(1)
        A1m[end, :] .= TR(0)
        A2m[end, :] .= TR(0)
        Anm[end, :] .= TR(0)
    end

    if hi_bc != BCType.Natural
        throw(ArgumentError("hi_bc not implemented yet"))
    end

    return QNMKerrRadialChebCache{TR}(;
        a=a, s=s, m=m, A0=A0m, A1=A1m, A2=A2m, An=Anm, Lm=Matrix{Complex{TR}}(undef, n, n)
    )
end

"""
    qnm_kerr_cheb_δ!(cache, ω, l; l_max=20)

Perform one iteration step in the QNM eigenvalue search.

# Arguments
- `cache`: Pre-computed QNMKerrRadialChebCache object
- `ω`: Current frequency guess
- `l`: Angular mode number
- `l_max`: Maximum l value for angular eigenvalue calculation (default: 20)

# Returns
- Difference between the seperation constant between the radial and angular equations
"""
function qnm_kerr_cheb_δ!(
    cache::QNMKerrRadialChebCache{TR}, ω::Complex{TR}, l::Integer; l_max::Integer=20
) where {TR<:AbstractFloat}
    @unpack_QNMKerrRadialChebCache cache

    c = a * ω
    ω2 = ω^2
    @.. Lm = A0 + ω * A1 + ω2 * A2
    Λ_vals = eigvals!(Lm, An; sortby=abs)

    A_idx = sws_eigvalidx(s, l, m)
    A = sws_eigvals(TR, s, c, m, l_max)[A_idx]

    Λ = -A
    δ = argmin(vi -> abs(vi - Λ), Λ_vals) - Λ

    return δ
end

export QNMKerrRadialChebCache, @unpack_QNMKerrRadialChebCache
export qnm_kerr_radial_cheb_cache, qnm_kerr_cheb_δ!
