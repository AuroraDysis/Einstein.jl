@kwdef struct QNMKerrContext{TR<:AbstractFloat,TI<:Integer}
    a::TR
    s::TI
    l::TI
    m::TI
    n::TI
    A_guess::Union{Complex{TR},Nothing}
    poles::Vector{Complex{TR}}
    l_max::TI
    cf_n_inv::TI
    cf_N_min::TI
    cf_N_max::TI
    cf_tol::TR
    angular_matrix::Matrix{Complex{TR}}
end

include("radial.jl")

function qnm_kerr_step!(
    ωvec::SVector{2,TR}, ctx::QNMKerrContext{TR,TI}
)::SVector{2,TR} where {TR<:AbstractFloat,TI<:Integer}
    (; a, s, l, m, l_max, A_guess, poles, angular_matrix) = ctx

    ω = ωvec[1] + im * ωvec[2]

    c = a * ω
    sws_eigM!(angular_matrix, s, c, m, l_max)
    A_vals = eigvals!(angular_matrix; sortby=abs)
    if isnothing(A_guess)
        A_idx = sws_eigvalidx(s, l, m)
        A = A_vals[A_idx]
    else
        A = argmin(vi -> abs(vi - A_guess), A_vals)
    end

    inv_cf, cf_error, cf_iter = qnm_kerr_radial(ctx, A, ω)

    pole_factors = prod(ω .- poles)
    supp_err = inv_cf / pole_factors

    return SA[real(supp_err), imag(supp_err)]
end

"""
    qnm_kerr(a, s, l, m, ω_guess; kwargs...)

Find the Kerr QNM using the Leaver's method for the radial equation and the Cook-Zalutskiy approach for the angular sector.

# Arguments
- `params`: QNMKerrCFParams object containing the Kerr parameters and initial guess
- `alg`: Nonlinear algorithm to use for the eigenvalue search (default: RobustMultiNewton with AutoFiniteDiff)
- `kwargs`: Additional keyword arguments to pass to the nonlinear solver
"""
function qnm_kerr(
    a::TF,
    s::TI,
    l::TI,
    m::TI,
    ω_guess::Complex{TF};
    n::TI=zero(TI),
    A_guess::Union{Complex{TF},Nothing}=nothing,
    poles::Vector{Complex{TF}}=Complex{TF}[],
    l_max::Union{TI,Nothing}=nothing,
    cf_n_inv::TI=TI(0),
    cf_N_min::TI=TI(300),
    cf_N_max::TI=TI(100000),
    cf_tol::TF=typetol(TF),
    nonlinear_algorithm::TA=RobustMultiNewton(; autodiff=AutoFiniteDiff()),
    angular_matrix_cache::Union{Matrix{Complex{TF}},Nothing}=nothing,
    kwargs...,
) where {TF<:AbstractFloat,TI<:Integer,TA<:AbstractNonlinearAlgorithm}
    if isnothing(l_max)
        l_max = l + 20
    end

    l_min = sws_l_min(s, m)
    l_size = l_max - l_min + 1
    angular_matrix = if isnothing(angular_matrix_cache)
        zeros(Complex{TF}, l_size, l_size)
    else
        angular_matrix_cache
    end

    ctx = QNMKerrContext{TF,TI}(;
        a=a,
        s=s,
        l=l,
        m=m,
        n=n,
        A_guess=A_guess,
        poles=poles,
        l_max=l_max,
        cf_n_inv=cf_n_inv,
        cf_N_min=cf_N_min,
        cf_N_max=cf_N_max,
        cf_tol=cf_tol,
        angular_matrix=angular_matrix,
    )

    ω0 = SVector(real(ω_guess), imag(ω_guess))
    prob = NonlinearProblem(qnm_kerr_step!, ω0, ctx)
    sol = solve(prob, nonlinear_algorithm; kwargs...)
    return Complex{TF}(sol.u...)
end

export qnm_kerr, qnm_kerr_step!, QNMKerrContext
