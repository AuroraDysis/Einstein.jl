@with_kw struct QNMKerrDIParams{TR<:AbstractFloat,TI<:Integer}
    a::TR
    s::TI
    l::TI
    m::TI
    n::TI
    ω_guess::Complex{TR}
    A_guess::Union{Complex{TR},Nothing} = nothing
    ρ_min::TR = TR(1//100)
    ρ_max::TR = one(TR) / (1 + sqrt(1 - a^2)) - TR(1//100)
    odealg::AbstractODEAlgorithm = Vern9()
    lo_bc::BCType.T = BCType.Natural
    hi_bc::BCType.T = BCType.Natural
    series_order::TI = 200
    abstol = typetol(TR)
    reltol = typetol(TR)
    l_max::TI = l + 20
end

struct QNMKerrDICache{TR<:AbstractFloat,TI<:Integer}
    params::QNMKerrDIParams{TR,TI}
    ω::Base.RefValue{Complex{TR}}
    Λ::Base.RefValue{Complex{TR}}
    M::Matrix{Complex{TR}}

    function QNMKerrDICache{TR,TI}(
        params::QNMKerrDIParams{TR,TI}
    ) where {TR<:AbstractFloat,TI<:Integer}
        @unpack_QNMKerrDIParams params

        ω = Ref{ComplexF64}()
        Λ = Ref{ComplexF64}()

        l_min = sws_l_min(s, m)
        l_size = l_max - l_min + 1
        M = zeros(Complex{TR}, l_size, l_size)

        return new{TR,TI}(params, ω, Λ, M)
    end
end

@doc raw"""
    qnm_kerr_radial_di_horizon_series(
        s::Integer,
        ρ::TR,
        a::TR,
        m::Integer,
        ω::Complex{TR},
        Λ::Complex{TR},
        series_order::Integer,
    ) where {TR<:AbstractFloat}

Calculate the series expansion of the radial equation at the horizon.

# Arguments
- `s::Integer`: spin weight of the field.
- `ρ::TR`: compactified radial coordinate.
- `a::TR`: black hole spin parameter.
- `m::Integer`: azimuthal mode number.
- `ω::Complex{TR}`: QNM frequency.
- `Λ::Complex{TR}`: separation constant.
- `series_order::Integer`: order of the series expansion.

# Returns
- `R0::Complex{TR}`: value of the radial function.
- `dR0::Complex{TR}`: value of the radial derivative.

# References
- [Ripley:2022ypi](@citet*)
"""
function qnm_kerr_radial_di_horizon_series(
    s::Integer,
    ρ::TR,
    a::TR,
    m::Integer,
    ω::Complex{TR},
    Λ::Complex{TR},
    series_order::Integer,
) where {TR<:AbstractFloat}
    # h[0] = 1, h[-1] (ρ - ρh) = 0, h[-2] (ρ - ρh)^2 = 0
    series_R_im1 = one(Complex{TR})
    series_R_im2 = zero(Complex{TR})
    series_R_im3 = zero(Complex{TR})

    R0 = one(Complex{TR})
    dR0 = zero(Complex{TR})

    ρh = one(TR) / (1 + sqrt(1 - a^2))

    ρmρh = ρ - ρh
    ρmρh2 = ρmρh * ρmρh
    ρmρh3 = ρmρh2 * ρmρh
    for i in 1:series_order
        coeff_im1 =
            (
                2 * a * m + 4im * i * (-1 + i + 1im * a * m) -
                1im * a^2 * (2 + i * (-5 + 5i - 2s) + 2s - Λ) * ρh +
                2im * (1 + s - Λ + 8im * ω) +
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

@doc raw"""
    qnm_kerr_radial_di_inf_series(
        s::Integer,
        ρ::TR,
        a::TR,
        m::Integer,
        ω::Complex{TR},
        Λ::Complex{TR},
        series_order::Integer,
    ) where {TR<:AbstractFloat}

Calculate the series expansion of the radial equation at the horizon.

# Arguments
- `s::Integer`: spin weight of the field.
- `ρ::TR`: compactified radial coordinate.
- `a::TR`: black hole spin parameter.
- `m::Integer`: azimuthal mode number.
- `ω::Complex{TR}`: QNM frequency.
- `Λ::Complex{TR}`: separation constant.
- `series_order::Integer`: order of the series expansion.

# Returns
- `R0::Complex{TR}`: value of the radial function.
- `dR0::Complex{TR}`: value of the radial derivative.

# References
- [Ripley:2022ypi](@citet*)
"""
function qnm_kerr_radial_di_inf_series(
    s::Integer,
    ρ::TR,
    a::TR,
    m::Integer,
    ω::Complex{TR},
    Λ::Complex{TR},
    series_order::Integer,
) where {TR<:AbstractFloat}
    # h[0] = 1, h[-1] ρ = 0, h[-2] ρ^2 = 0
    series_R_im1 = one(Complex{TR})
    series_R_im2 = zero(Complex{TR})
    series_R_im3 = zero(Complex{TR})

    R0 = one(Complex{TR})
    dR0 = zero(Complex{TR})

    ρ2 = ρ * ρ
    ρ3 = ρ2 * ρ
    for i in 1:series_order
        coeff_im1 =
            im * (
                2 * s - i * (-1 + i + 2s) - Λ +
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

@doc raw"""
    qnm_kerr_radial_di_rhs(
        u::SVector{2,Complex{TR}},
        cache::QNMKerrDICache{TR,TI},
        ρ::TR,
    ) where {TR<:AbstractFloat,TI<:Integer}

Calculate the right-hand side of the radial equation in hyperboloidal coordinates.

# Arguments
- `u::SVector{2,Complex{TR}}`: radial function and its derivative.
- `cache::QNMKerrDICache{TR,TI}`: cache object.
- `ρ::TR`: compactified radial coordinate.

# References
- [Ripley:2022ypi](@citet*)
"""
function qnm_kerr_radial_di_rhs(
    u::SVector{2,Complex{TR}}, cache::QNMKerrDICache{TR,TI}, ρ::TR
) where {TR<:AbstractFloat,TI<:Integer}
    @unpack_QNMKerrDIParams cache.params

    ω = cache.ω[]
    Λ = cache.Λ[]

    R, dR = u
    ddR =
        (
            R * (
                -Λ +
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

@doc raw"""
    qnm_kerr_radial_di_δ!(
        x::SVector{2,TR},
        cache::QNMKerrDICache{TR,TI},
    ) where {TR<:AbstractFloat,TI<:Integer}

Integrate the radial equation from infinity (ρ_min) to the horizon (ρ_max) and calculate the difference compared to
- series expansion at the horizon (Natural boundary condition)
- 0 (Dirichlet boundary condition)

# Arguments
- `x::SVector{2,TR}`: real and imaginary parts of the QNM frequency.
- `cache::QNMKerrDICache{TR,TI}`: cache object.

# Returns
- `δ::SVector{2,TR}`: real and imaginary parts of the difference.

# References
- [Ripley:2022ypi](@citet*)
"""
function qnm_kerr_radial_di_δ!(
    cache::QNMKerrDICache{TR,TI}
)::Complex{TR} where {TR<:AbstractFloat,TI<:Integer}
    @unpack_QNMKerrDIParams cache.params

    ω = cache.ω[]
    Λ = cache.Λ[]

    u0 = qnm_kerr_radial_di_inf_series(s, ρ_min, a, m, ω, Λ, series_order)
    tspan = (ρ_min, ρ_max)
    prob = ODEProblem(qnm_kerr_radial_di_rhs, u0, tspan, cache)

    sol = solve(
        prob,
        odealg;
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
            Rend = sol[end][1]
            return Rend
        else
            sol_end = sol.u[end]
            sol_end_match = qnm_kerr_radial_di_horizon_series(
                s, ρ_max, a, m, ω, Λ, series_order
            )
            δ = sol_end[2] ./ sol_end[1] .- sol_end_match[2] ./ sol_end_match[1]
            return δ
        end

        if hi_bc != BCType.Natural
            throw(ArgumentError("hi_bc not implemented yet"))
        end
    end
end

function qnm_kerr_di_δ!(
    ωvec::SVector{2,TR}, cache::QNMKerrDICache{TR,TI}
)::SVector{2,TR} where {TR<:AbstractFloat,TI<:Integer}
    @unpack_QNMKerrDIParams cache.params

    ω = Complex{TR}(ωvec[1], ωvec[2])

    M = cache.M
    c = a * ω
    sws_eigM!(M, s, c, m, l_max)
    A_vals = eigvals!(M; sortby=abs)
    if isnothing(A_guess)
        A_idx = sws_eigvalidx(s, l, m)
        A = A_vals[A_idx]
    else
        A = argmin(vi -> abs(vi - A_guess), A_vals)
    end

    Λ = -A

    cache.ω[] = ω
    cache.Λ[] = Λ

    δ = qnm_kerr_radial_di_δ!(cache)

    return SA[real(δ), imag(δ)]
end

"""
    qnm_kerr_di(params::QNMKerrDIParams{TR}, alg::AbstractNonlinearAlgorithm=RobustMultiNewton(autodiff=AutoFiniteDiff()); kwargs...) where {TR<:AbstractFloat}

Find the Kerr QNM using the direct integration method for the radial equation and the Cook-Zalutskiy approach for the angular sector.

# Arguments
- `params`: QNMKerrDIParams  object containing the Kerr parameters and initial guess
- `alg`: Nonlinear algorithm to use for the eigenvalue search (default: RobustMultiNewton with AutoFiniteDiff)
- `kwargs`: Additional keyword arguments to pass to the nonlinear solver
"""
function qnm_kerr_di(
    params::QNMKerrDIParams{TR,TI},
    alg::Union{AbstractNonlinearAlgorithm,Nothing}=RobustMultiNewton(;
        autodiff=AutoFiniteDiff()
    );
    kwargs...,
) where {TR<:AbstractFloat,TI<:Integer}
    @unpack_QNMKerrDIParams params

    cache = QNMKerrDICache{TR,TI}(params)
    ω0 = SA[real(ω_guess), imag(ω_guess)]
    prob = NonlinearProblem(qnm_kerr_di_δ!, ω0, cache)
    sol = solve(prob, alg; kwargs...)
    ω = Complex{TR}(sol.u[1], sol.u[2])

    return ω
end

export qnm_kerr_radial_di_horizon_series,
    qnm_kerr_radial_di_inf_series,
    qnm_kerr_radial_di_rhs,
    QNMKerrDICache,
    QNMKerrDIParams,
    @unpack_QNMKerrDIParams,
    qnm_kerr_radial_di_δ!,
    qnm_kerr_di
