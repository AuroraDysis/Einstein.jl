@with_kw struct QNMKerrCache{TR<:AbstractFloat}
    a::TR # black hole spin
    s::Integer
    m::Integer
    A0::AbstractMatrix{Complex{TR}}
    A1::AbstractMatrix{Complex{TR}}
    A2::AbstractMatrix{TR}
    An::AbstractMatrix{Complex{TR}}
    has_v::Bool
    v::Vector{Complex{TR}}
    Lm::Matrix{Complex{TR}}
end

function qnm_kerrpep_cache(
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

    dom = zero(TR) .. one(TR)
    chebSpace = Chebyshev(dom)
    ultraSpace = Ultraspherical(2, dom)
    conversion = Conversion(chebSpace, ultraSpace)
    conversionA1 = Conversion(Ultraspherical(1, dom), ultraSpace)

    ρn = Fun(chebSpace)
    ρ = ρn * (ρ_max - ρ_min) + ρ_min
    c02 = -ρ^2 * (1 - 2 * M * ρ + a^2 * ρ^2) / (ρ_max - ρ_min)^2
    c01 =
        -2 * ρ * (1 + s - (1im * a * m + M * (3 + s)) * ρ + 2 * a^2 * ρ^2) / (ρ_max - ρ_min)
    c00 = 2 * ρ * (1im * a * m + M * (1 + s) - a^2 * ρ)
    A0c = (c02 * 𝒟^2 + c01 * 𝒟 + c00):chebSpace

    c11 = 2im * (1 + ρ^2 * (-8 * M^2 + a^2 * (1 + 4 * M * ρ))) / (ρ_max - ρ_min)
    c10 =
        2im * a^2 * ρ * (1 + 6 * M * ρ) +
        2 * a * m * (1 + 4 * M * ρ) +
        4im * M * (s - 2 * M * (s + 2) * ρ)
    A1 = (c11 * 𝒟 + c10):chebSpace
    A1c = conversionA1 * A1

    c20 = -16 * M^2 * (1 + 2 * M * ρ) + a^2 * (1 + 4 * M * ρ)^2
    A2 = (c20 * conversion):chebSpace

    A0m = A0c[1:n, 1:n]
    A1m = A1c[1:n, 1:n]
    A2m = A2[1:n, 1:n]
    Anm = Matrix{Complex{TR}}(conversion[1:n, 1:n])

    if lo_bc == BCType.Dirichlet
        A0m[end, :] .= TR(1)
        A1m[end, :] .= TR(0)
        A2m[end, :] .= TR(0)
        Anm[end, :] .= TR(0)
    end

    if hi_bc != BCType.Natural
        throw(ArgumentError("hi_bc not implemented yet"))
    end

    return QNMKerrCache{TR}(;
        a=a,
        s=s,
        m=m,
        A0=A0m,
        A1=A1m,
        A2=A2m,
        An=Anm,
        has_v=false,
        v=Vector{Complex{TR}}(undef, n),
        Lm=Matrix{Complex{TR}}(undef, n, n),
    )
end

function qnm_kerrpep_step!(
    cache::QNMKerrCache{TR}, ω::Complex{TR}
) where {TR<:AbstractFloat}
    @unpack_QNMKerrCache cache

    @.. Lm = A0 + ω * A1 + ω^2 * A2
    return eigvals!(Lm, An; sortby=abs)
end

export QNMKerrCache, @unpack_QNMKerrCache
export qnm_kerrpep_cache
