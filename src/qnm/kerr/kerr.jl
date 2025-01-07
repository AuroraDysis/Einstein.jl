include("radial.jl")
include("cheb.jl")
include("di.jl")

@with_kw struct QNMKerrParams{TR<:AbstractFloat}
    a::TR
    s::Integer
    l::Integer
    m::Integer
    n::Integer = 0
    ω_guess::Complex{TR}
    A_guess::Complex{TR} = sws_A0(s, l)
    poles::Vector{Complex{TR}} = []
    l_max::Integer = l + 20
    cf_N_min::Integer = 300
    cf_N_max::Integer = 100000
    cf_tol::TR = typetol(TR)
end

struct QNMKerrCache{TR<:AbstractFloat}
    params::QNMKerrParams{TR}
    M::Matrix{Complex{TR}}

    function QNMKerrCache{TR}(params::QNMKerrParams{TR}) where {TR<:AbstractFloat}
        @unpack_QNMKerrParams params

        l_min = sws_l_min(s, m)
        l_size = l_max - l_min + 1
        M = zeros(Complex{TR}, l_size, l_size)
        return new{TR}(params, M)
    end
end

function qnm_kerrnep_δ(
    x::SVector{2,TR}, cache::QNMKerrCache{TR}
)::SVector{2,TR} where {TR<:AbstractFloat}
    @unpack_QNMKerrParams cache.params
    M = cache.M

    ω = x[1] + im * x[2]
    c = a * ω

    sws_eigM!(M, s, c, m, l_max)
    evals = eigvals!(M; sortby=abs)
    evals_idx = sws_eigvalidx(s, l, m)
    A = evals[evals_idx]

    inv_cf, cf_error, cf_iter = qnm_kerr_radial_cf(TR, a, s, m, A, ω)

    pole_factors = prod(ω .- poles)
    supp_err = inv_cf / pole_factors

    return SA[real(supp_err), imag(supp_err)]
end

"""
    qnm_kerr(params::QNMKerrParams{TR}; alg=RobustMultiNewton(autodiff=AutoFiniteDiff()), kwargs...)

Find the Kerr QNM using the Leaver's method for the radial equation and the Cook-Zalutskiy approach for the angular sector.

# Arguments
- `params`: QNMKerrParams object containing the Kerr parameters and initial guess
- `alg`: Nonlinear algorithm to use for the eigenvalue search (default: RobustMultiNewton with AutoFiniteDiff)
- `kwargs`: Additional keyword arguments to pass to the nonlinear solver
"""
function qnm_kerr(
    params::QNMKerrParams{TR},
    alg::Union{AbstractNonlinearAlgorithm,Nothing}=RobustMultiNewton(;
        autodiff=AutoFiniteDiff()
    ),
    kwargs...,
) where {TR<:AbstractFloat}
    @unpack_QNMKerrParams params

    cache = QNMKerrCache{TR}(params)
    ω0 = SA[real(ω_guess), imag(ω_guess)]
    prob = NonlinearProblem(qnm_kerrnep_δ, ω0, cache)

    sol = solve(prob, alg, kwargs...)
    ω = sol.u[1] + im * sol.u[2]
    return ω
end

export QNMKerrParams, QNMKerrCache, qnm_kerrnep_δ, qnm_kerr, @unpack_QNMKerrParams
