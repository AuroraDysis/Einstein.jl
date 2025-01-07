function qnm_kerr_radial_di_horizon_series(
    s::Integer, ρ::T, a::T, m::Integer, ω::Complex{T}, A::Complex{T}, series_order::Integer
) where {T<:AbstractFloat}
    # h[0] = 1, h[-1] (ρ - ρh) = 0, h[-2] (ρ - ρh)^2 = 0
    series_R_im1 = one(Complex{T})
    series_R_im2 = zero(Complex{T})
    series_R_im3 = zero(Complex{T})

    R0 = one(Complex{T})
    dR0 = zero(Complex{T})

    ρh = one(T) / (1 + sqrt(1 - a^2))

    ρmρh = ρ - ρh
    ρmρh2 = ρmρh * ρmρh
    ρmρh3 = ρmρh2 * ρmρh
    for i in 1:series_order
        coeff_im1 =
            (
                2 * a * m + 4im * i * (-1 + i + 1im * a * m) -
                1im * a^2 * (2 + i * (-5 + 5i - 2s) + 2s - A) * ρh +
                2im * (1 + s - A + 8im * ω) +
                2 *
                (16i + a * (a - 2a * i + 6im * m + a * (6 - 12i - im * a * m + 2s) * ρh)) *
                ω - im * (64 + a^4 * ρh - 2a^2 * (5 + 16ρh)) * ω^2
            ) / ρh^2

        coeff_im2 =
            (
                2 *
                (-1 + i - 4im * ω) *
                (
                    2im * (-1 + i + s) +
                    8ω +
                    a * (
                        -2m + a^2 * m * ρh - 1im * a * (-2 + (-1 + s) * ρh + i * (2 + ρh)) +
                        a^3 * ρh * ω - 2a * (3 + 2ρh) * ω
                    )
                )
            ) / ρh^2

        coeff_im3 =
            (a^2 * (-2 + a^2 * ρh) * (-1 + i - 4im * ω) * (1im * (-2 + i) + 4ω)) / ρh^2

        den =
            (
                2 *
                i *
                (a * m + 1im * s - 1im * (i + a^2 * (-i + s) * ρh) - 4ω + 2 * a^2 * ρh * ω)
            ) / ρh

        series_R_i =
            (
                ρmρh * coeff_im1 * series_R_im1 +
                ρmρh2 * coeff_im2 * series_R_im2 +
                ρmρh3 * coeff_im3 * series_R_im3
            ) / den

        R0 += series_R_i
        dR0 += i * series_R_i / ρmρh

        series_R_im3 = series_R_im2
        series_R_im2 = series_R_im1
        series_R_im1 = series_R_i
    end

    return SA[R0, dR0]
end

function qnm_kerr_radial_di_inf_series(
    s::Integer, ρ::T, a::T, m::Integer, ω::Complex{T}, A::Complex{T}, series_order::Integer
) where {T<:AbstractFloat}
    # h[0] = 1, h[-1] ρ = 0, h[-2] ρ^2 = 0
    series_R_im1 = one(Complex{T})
    series_R_im2 = zero(Complex{T})
    series_R_im3 = zero(Complex{T})

    R0 = one(Complex{T})
    dR0 = zero(Complex{T})

    ρ2 = ρ * ρ
    ρ3 = ρ2 * ρ
    for i in 1:series_order
        coeff_im1 =
            im * (
                2 * s - i * (-1 + i + 2s) - A +
                2 * a * m * ω +
                4im * s * ω +
                (-16 + a^2) * ω^2
            )
        coeff_im2 =
            2 * (1im * (i - 1) + 4 * ω) * (-1 + i + s - 4im * ω + im * a * (m + a * ω))
        coeff_im3 = -im * a^2 * (-2 + i - 4im * ω) * (-1 + i - 4im * ω)
        den = 2 * i * ω
        series_R_i =
            (
                ρ * coeff_im1 * series_R_im1 +
                ρ2 * coeff_im2 * series_R_im2 +
                ρ3 * coeff_im3 * series_R_im3
            ) / den

        R0 += series_R_i
        dR0 += i * series_R_i / ρ

        series_R_im3 = series_R_im2
        series_R_im2 = series_R_im1
        series_R_im1 = series_R_i
    end

    return SA[R0, dR0]
end

@with_kw struct QNMKerrRadialDIParams{T<:AbstractFloat}
    a::T
    s::Integer
    m::Integer
    ω_guess::Complex{T}
    A_guess::Complex{T} = sws_A0(s, m)
    ρ_min::T = T(1//100)
    ρ_max::T = one(T) / (1 + sqrt(1 - a^2)) - T(1//100)
    alg::AbstractODEAlgorithm = Vern9()
    lo_bc::BCType.T = BCType.Natural
    hi_bc::BCType.T = BCType.Natural
    series_order::Integer = 200
    abstol = typetol(T)
    reltol = typetol(T)
end

mutable struct QNMKerrRadialDICache{T<:AbstractFloat}
    params::QNMKerrRadialDIParams{T}
    ω::Complex{T}
    A::Complex{T}
end

function qnm_kerr_radial_di_rhs(
    u::SVector{2,Complex{T}}, cache::QNMKerrRadialDICache{T}, ρ::T
) where {T<:AbstractFloat}
    @unpack_QNMKerrRadialDIParams cache.params

    ω = cache.ω
    A = cache.A

    R, dR = u
    ddR =
        (
            R * (
                -A +
                ω * (2a * m + 4im * s - 16ω + a^2 * ω) +
                2a^2 * ρ^2 * (-1 + 6im * ω + 8ω^2) +
                2ρ * (im + 4ω) * (-im * (1 + s) - 4ω + a * (m + a * ω))
            ) +
            dR * (
                2im * ω +
                2ρ * (
                    -1 +
                    s * (-1 + ρ) +
                    ρ * (3 - 8im * ω + im * a * (m + a * (2im * ρ + ω + 4ρ * ω)))
                )
            )
        ) / (ρ^2 * (ρ * (a^2 * ρ - 2) + 1))
    return SA[dR, ddR]
end

function qnm_kerr_radial_di_δ(
    x::SVector{2,TR}, cache::QNMKerrRadialDICache{TR}
)::SVector{2,Complex{TR}} where {TR<:AbstractFloat}
    @unpack_QNMKerrRadialDIParams cache.params

    ω = Complex{TR}(x[1], x[2])
    cache.ω = ω
    A = cache.A

    u0 = qnm_kerr_radial_di_inf_series(s, ρ_min, a, m, ω, A, series_order)
    tspan = (ρ_min, ρ_max)
    prob = ODEProblem(qnm_kerr_radial_di_rhs, u0, tspan, cache)

    sol = solve(
        prob,
        alg;
        abstol=abstol,
        reltol=reltol,
        save_everystep=false,
        save_start=false,
        save_end=true,
    )

    if sol.retcode != :Success
        error("ODE solver failed")
    else
        if lo_bc == BCType.Dirichlet
            sol_end = sol[end]
            return sol_end
        else
            sol_end = sol.u[end]
            sol_end_match = qnm_kerr_radial_di_horizon_series(s, ρ_max, a, m, ω, A, series_order)
            δ = sol_end ./ sol_end[1] .- sol_end_match ./ sol_end_match[1]
            return δ
        end

        if hi_bc != BCType.Natural
            throw(ArgumentError("hi_bc not implemented yet"))
        end
    end
end

export qnm_kerr_radial_di_horizon_series,
    qnm_kerr_radial_di_inf_series,
    qnm_kerr_radial_di_rhs,
    QNMKerrRadialDICache,
    QNMKerrRadialDIParams,
    @unpack_QNMKerrRadialDIParams,
    qnm_kerr_radial_di_δ
