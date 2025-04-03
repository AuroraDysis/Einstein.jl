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

@doc raw"""
    qnm_kerr(a, s, l, m, ω_guess; kwargs...)

Find the Kerr QNM using Leaver's method for the radial equation and the Cook-Zalutskiy approach for the angular sector.
We use the unit $M = G = c = 1$.

# Arguments
- `a`: Kerr parameter (dimensionless spin parameter, 0 ≤ |a| ≤ 1)
- `s`: Spin weight of the field (±2 for gravitational perturbations, ±1 for electromagnetic perturbations, 0 for scalar perturbations)
- `l`: Angular momentum quantum number (l ≥ |s|)
- `m`: Azimuthal harmonic index (-l ≤ m ≤ l)
- `ω_guess`: Initial guess for the QNM frequency

# Keyword Arguments
- `n`: Overtone number (default: 0 for fundamental mode)
- `A_guess`: Initial guess for the angular separation constant (default: nothing)
- `poles`: Vector of known poles to subtract from the continued fraction (default: empty)
- `l_max`: Maximum angular momentum for angular eigenvalue calculation (default: l + 20)
- `cf_n_inv`: Number of inverse terms in the continued fraction (default: 0)
- `cf_N_min`: Minimum number of terms in the continued fraction (default: 300)
- `cf_N_max`: Maximum number of terms in the continued fraction (default: 100000)
- `cf_tol`: Tolerance for continued fraction convergence (default: machine epsilon)
- `nonlinear_algorithm`: Algorithm for solving the nonlinear eigenvalue problem (default: RobustMultiNewton with AutoFiniteDiff)
- `angular_matrix_cache`: Pre-allocated matrix for angular eigenvalue calculation (default: nothing)
- `kwargs...`: Additional keyword arguments passed to the nonlinear solver

# Returns
- `Complex{TF}`: The complex QNM frequency $\omega$

# Examples
```julia
a = 0.7
s = 2
l = 2
m = 2
n = 0
ω = 0.532600243551018 - 0.08079287315500766im
l_max = 20

ω_pert = ω + rand(Complex{Float64}) / 1000

ωsol = qnm_kerr(a, s, l, m, ω_pert; n=n, l_max=l_max)
```

# Notes
- The function uses a combination of Leaver's method for the radial equation and the Cook-Zalutskiy approach for the angular sector
- The continued fraction method is used to solve the radial equation
- The angular eigenvalue problem is solved using matrix eigenvalue methods
- The function returns the complex frequency where the imaginary part is negative (damped modes)
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
