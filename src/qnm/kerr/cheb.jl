@with_kw struct QNMKerrChebParams{TR<:AbstractFloat}
    a::TR
    s::Integer
    l::Integer
    m::Integer
    n::Integer
    cheb_n::Integer
    ω_guess::Complex{TR}
    A_guess::Complex{TR} = sws_A0(s, l)
    l_max::Integer = l + 20
    ρ_min::TR = zero(TR)
    ρ_max::TR = one(TR)
    lo_bc::BCType.T = BCType.Natural
    hi_bc::BCType.T = BCType.Natural
end

struct QNMKerrChebCache{TR<:AbstractFloat}
    params::QNMKerrChebParams{TR}
    A0::AbstractMatrix{Complex{TR}}
    A1::AbstractMatrix{Complex{TR}}
    A2::AbstractMatrix{TR}
    An::AbstractMatrix{Complex{TR}}
    LA::AbstractMatrix{Complex{TR}}
    LB::AbstractMatrix{Complex{TR}}
    M::AbstractMatrix{Complex{TR}}

    function QNMKerrChebCache{TR}(params::QNMKerrChebParams{TR}) where {TR<:AbstractFloat}
        @unpack_QNMKerrChebParams params

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

        A0m = Matrix(@view(A0c[1:cheb_n, 1:cheb_n]))
        A1m = Matrix(@view(A1c[1:cheb_n, 1:cheb_n]))
        A2m = Matrix(@view(A2[1:cheb_n, 1:cheb_n]))
        Anm = Matrix(@view(conversion[1:cheb_n, 1:cheb_n]))

        if lo_bc == BCType.Dirichlet
            A0m[end, :] .= TR(1)
            A1m[end, :] .= TR(0)
            A2m[end, :] .= TR(0)
            Anm[end, :] .= TR(0)
        end

        if hi_bc != BCType.Natural
            throw(ArgumentError("hi_bc not implemented yet"))
        end

        l_min = sws_l_min(s, m)
        l_size = l_max - l_min + 1
        M = zeros(Complex{TR}, l_size, l_size)

        LA = zeros(Complex{TR}, cheb_n, cheb_n)
        LB = zeros(Complex{TR}, cheb_n, cheb_n)

        return new{TR}(params, A0m, A1m, A2m, Anm, LA, LB, M)
    end
end

"""
    qnm_kerr_cheb_δ!(cache, ω, l; l_max=20)

Perform one iteration step in the QNM eigenvalue search.

# Arguments
- `cache`: Pre-computed QNMKerrChebCache object
- `ω`: Current frequency guess
- `l`: Angular mode number
- `l_max`: Maximum l value for angular eigenvalue calculation (default: 20)

# Returns
- Difference between the seperation constant between the radial and angular equations
"""
function qnm_kerr_cheb_δ!(
    ωvec::SVector{2,TR}, cache::QNMKerrChebCache{TR}
)::SVector{2,TR} where {TR<:AbstractFloat}
    @unpack params, A0, A1, A2, An, LA, LB, M = cache
    @unpack_QNMKerrChebParams params

    ω = Complex{TR}(ωvec[1], ωvec[2])

    ω2 = ω^2
    @.. LA = A0 + ω * A1 + ω2 * A2
    LB .= An
    Λ_vals = eigvals!(LA, LB; sortby=abs)

    c = a * ω
    sws_eigM!(M, s, c, m, l_max)
    A_vals = eigvals!(M; sortby=abs)
    A_idx = sws_eigvalidx(s, l, m)
    A = A_vals[A_idx]

    Λ = -A
    δ = argmin(vi -> abs(vi - Λ), Λ_vals) - Λ

    return SA[real(δ), imag(δ)]
end

function qnm_kerr_cheb(
    params::QNMKerrChebParams{TR},
    alg::Union{AbstractNonlinearAlgorithm,Nothing}=RobustMultiNewton(;
        autodiff=AutoFiniteDiff()
    );
    kwargs...,
) where {TR<:AbstractFloat}
    @unpack_QNMKerrChebParams params

    cache = QNMKerrChebCache{TR}(params)
    ω0 = SA[real(ω_guess), imag(ω_guess)]
    prob = NonlinearProblem(qnm_kerr_cheb_δ!, ω0, cache)
    sol = solve(prob, alg; kwargs...)
    ω = Complex{TR}(sol.u[1], sol.u[2])

    return ω
end

export QNMKerrChebParams, QNMKerrChebCache, @unpack_QNMKerrChebParams
export qnm_kerr_cheb_δ!, qnm_kerr_cheb
