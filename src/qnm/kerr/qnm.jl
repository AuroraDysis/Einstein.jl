@kwdef struct QNMKerrContext{TR<:AbstractFloat,TI<:Integer}
    a::TR
    s::TI
    l::TI
    m::TI
    n::TI = 0
    poles::Vector{Complex{TR}}
    l_max::TI
    cf_N_min::TI
    cf_N_max::TI
    cf_tol::TR
end

include("radial.jl")

function qnm_kerr(
    a::TF,
    s::TI,
    l::TI,
    m::TI;
    n::TI=zero(TI),
    ω_guess::Complex{TF},
    A_guess::Union{Complex{TF},Nothing}=nothing,
    poles::Vector{Complex{TF}}=Complex{TF}[],
    l_max::TI=l + 20,
    cf_N_min::TI=300,
    cf_N_max::TI=100000,
    cf_tol::TF=typetol(TF),
) where {TF<:AbstractFloat,TI<:Integer}
    return ctx = QNMKerrContext{TF,TI}(;
        a=a,
        s=s,
        l=l,
        m=m,
        n=n,
        ω_guess=ω_guess,
        A_guess=A_guess,
        poles=poles,
        l_max=l_max,
        cf_N_min=cf_N_min,
        cf_N_max=cf_N_max,
        cf_tol=cf_tol,
    )
end

struct QNMKerrCFCache{TR<:AbstractFloat,TI<:Integer}
    params::QNMKerrCFParams{TR,TI}
    M::Matrix{Complex{TR}}

    function QNMKerrCFCache{TR,TI}(
        params::QNMKerrCFParams{TR,TI}
    ) where {TR<:AbstractFloat,TI<:Integer}
        @unpack_QNMKerrCFParams params

        l_min = sws_l_min(s, m)
        l_size = l_max - l_min + 1
        M = zeros(Complex{TR}, l_size, l_size)
        return new{TR,TI}(params, M)
    end
end



export qnm_kerr
