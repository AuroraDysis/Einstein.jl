include("radial.jl")
include("cheb.jl")

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
    cf_tol::TR = typetol(typeof(a))
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

function qnm_kerrnepδ(
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

    inv_cf, cf_error, cf_iter = qnm_kerrrad(TR, a, s, m, A, ω)

    pole_factors = prod(ω .- poles)
    supp_err = inv_cf / pole_factors

    return SA[real(supp_err), imag(supp_err)]
end

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
    prob = NonlinearProblem(qnm_kerrnepδ, ω0, cache)

    sol = solve(prob, alg, kwargs...)
    ω = sol.u[1] + im * sol.u[2]
    return ω
end

export QNMKerrParams, QNMKerrCache, qnm_kerrnepδ, qnm_kerr
