@with_kw struct QNMKerrCFParams{TR<:AbstractFloat}
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

struct QNMKerrCFCache{TR<:AbstractFloat}
    params::QNMKerrCFParams{TR}
    M::Matrix{Complex{TR}}

    function QNMKerrCFCache{TR}(params::QNMKerrCFParams{TR}) where {TR<:AbstractFloat}
        @unpack_QNMKerrCFParams params

        l_min = sws_l_min(s, m)
        l_size = l_max - l_min + 1
        M = zeros(Complex{TR}, l_size, l_size)
        return new{TR}(params, M)
    end
end

"""
MIT License

Copyright (c) 2025 Zhen Zhong
Copyright (c) 2019 Leo C. Stein

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

@doc raw"""
    qnm_kerr_radial_cf(
        ::Type{TR},
        a::TR,
        s::Integer,
        m::Integer,
        A::Complex{TR},
        ω::Complex{TR};
        n_inv::Integer=0,
        cf_tol::TR=typetol(TR),
        cf_N_min::Integer=300,
        cf_N_max::Integer=100000,
    ) where {TR<:AbstractFloat}

Calculate the radial function using the Leaver scheme [Cook:2014cta](@cite).

## Arguments
- `TR`: Type of the floating-point number.
- a::TR: Black hole spin.
- s::Integer: spin weight of the field.
- m::Integer: Azimuthal number.
- A::Complex{TR}: Angular separation constant.
- ω::Complex{TR}: QNM frequency.
- n_inv::Integer: Number of inversions.
- cf_tol::TR: Tolerance for convergence.
- cf_N_min::Integer: Minimum number of iterations.
- cf_N_max::Integer: Maximum number of iterations.

## Returns

- `inv_cf::Complex{TR}`: The inverse continued fraction.
- `error::TR`: The error.
- `iter::Integer`: The number of iterations.

## Leaver scheme

```math
\begin{aligned}
& D_0=\delta=1+s+2 \xi \\
& D_1=4 p-2 \alpha+\gamma-\delta-2 \\
& D_2=2 \alpha-\gamma+2 \\
& D_3=\alpha(4 p-\delta)-\sigma, \\
& D_4=\alpha(\alpha-\gamma+1)
\end{aligned}
```

```math
\begin{aligned}
\alpha_n & \equiv n^2+\left(D_0+1\right) n+D_0 \\
\beta_n & \equiv-2 n^2+\left(D_1+2\right) n+D_3 \\
\gamma_n & \equiv n^2+\left(D_2-3\right) n+D_4-D_2+2
\end{aligned}
```

The general continued fraction:
```math
0=\beta_0-\frac{\alpha_0 \gamma_1}{\beta_1-} \frac{\alpha_1 \gamma_2}{\beta_2-} \frac{\alpha_2 \gamma_3}{\beta_3-} \ldots
```

The truncated version of the continued fraction:
```math
\operatorname{Cf}(\mathrm{N}) \equiv \beta_0-\frac{\alpha_0 \gamma_1}{\beta_1-} \frac{\alpha_1 \gamma_2}{\beta_2-} \frac{\alpha_2 \gamma_3}{\beta_3-} \ldots \frac{\alpha_{\mathrm{N}-1} \gamma_{\mathrm{N}}}{\beta_{\mathrm{N}}+\alpha_{\mathrm{N}} \mathrm{r}_{\mathrm{N}}}
```

The nth inversion of this truncated continued fraction:
```math
\begin{aligned}
\operatorname{Cf}(\mathrm{n} ; \mathrm{N}) \equiv & \beta_n-\frac{\alpha_{n-1} \gamma_n}{\beta_{n-1}-} \frac{\alpha_{n-2} \gamma_{n-1}}{\beta_{n-2}-} \ldots \frac{\alpha_0 \gamma_1}{\beta_0} \\
& -\frac{\alpha_n \gamma_{n+1}}{\beta_{n+1}-} \frac{\alpha_{n+1} \gamma_{n+2}}{\beta_{n+2}-} \cdots \frac{\alpha_{N-1} \gamma_N}{\beta_N+\alpha_N r_N}
\end{aligned}
```
with $\operatorname{Cf}(0 ; \mathrm{N}) \equiv \operatorname{Cf}(\mathrm{N})$ and $0 \leq n < N$

## References
- [Cook:2014cta](@citet*)
- [qnm/qnm/radial.py at master · duetosymmetry/qnm](https://github.com/duetosymmetry/qnm/blob/master/qnm/radial.py)
"""
function qnm_kerr_radial_cf(
    ::Type{TR},
    a::TR,
    s::Integer,
    m::Integer,
    A::Complex{TR},
    ω::Complex{TR};
    n_inv::Integer=0,
    cf_tol::TR=typetol(TR),
    cf_N_min::Integer=300,
    cf_N_max::Integer=100000,
) where {TR<:AbstractFloat}
    root = sqrt(1 - a * a)
    r₊ = 1 + root
    r₋ = 1 - root

    σ₊ = (2 * ω * r₊ - m * a) / (2 * root)
    σ₋ = (2 * ω * r₋ - m * a) / (2 * root)

    ζ = 1im * ω # ζ₊
    ξ = -s - 1im * σ₊ # ξ₋
    η = -1im * σ₋ # η₊

    p = root * ζ
    α = 1 + s + ξ + η - 2 * ζ + s
    γ = 1 + s + 2 * η
    δ = 1 + s + 2 * ξ
    σ =
        A + a * a * ω * ω - 8 * ω * ω +
        p * (2 * α + γ - δ) +
        (1 + s - (γ + δ) / 2) * (s + (γ + δ) / 2)

    D0 = δ
    D1 = 4 * p - 2 * α + γ - δ - 2
    D2 = 2 * α - γ + 2
    D3 = α * (4 * p - δ) - σ
    D4 = α * (α - γ + 1)

    n = 0:n_inv
    αn = zeros(Complex{TR}, length(n))
    βn = zeros(Complex{TR}, length(n))
    γn = zeros(Complex{TR}, length(n))

    @.. αn = n * n + (D0 + 1) * n + D0
    @.. βn = -2 * n * n + (D1 + 2) * n + D3
    @.. γn = n * n + (D2 - 3) * n + D4 - D2 + 2

    conv1 = zero(Complex{TR})
    @inbounds for i in 1:n_inv
        conv1 = αn[i] / (βn[i] - γn[i] * conv1)
    end

    function rad_a(i)
        n = i + n_inv - 1
        return -(n * n + (D0 + 1) * n + D0) / (n * n + (D2 - 3) * n + D4 - D2 + 2)
    end

    function rad_b(i)
        if i == 0
            return 0
        end
        n = i + n_inv
        return (-2 * n * n + (D1 + 2) * n + D3) / (n * n + (D2 - 3) * n + D4 - D2 + 2)
    end

    conv2, error, iter = contfrac_lentz(TR, rad_a, rad_b, cf_tol, cf_N_min, cf_N_max)

    inv_cf = βn[n_inv + 1] - γn[n_inv + 1] * conv1 + γn[n_inv + 1] * conv2

    return inv_cf, error, iter
end

function qnm_kerr_cf_δ!(
    ωvec::SVector{2,TR}, cache::QNMKerrCFCache{TR}
)::SVector{2,TR} where {TR<:AbstractFloat}
    @unpack_QNMKerrCFParams cache.params
    M = cache.M

    ω = ωvec[1] + im * ωvec[2]
    c = a * ω

    sws_eigM!(M, s, c, m, l_max)
    A_vals = eigvals!(M; sortby=abs)
    A_idx = sws_eigvalidx(s, l, m)
    A = A_vals[A_idx]

    inv_cf, cf_error, cf_iter = qnm_kerr_radial_cf(TR, a, s, m, A, ω)

    pole_factors = prod(ω .- poles)
    supp_err = inv_cf / pole_factors

    return SA[real(supp_err), imag(supp_err)]
end

"""
    qnm_kerr_cf(params::QNMKerrCFParams{TR}; alg=RobustMultiNewton(autodiff=AutoFiniteDiff()), kwargs...)

Find the Kerr QNM using the Leaver's method for the radial equation and the Cook-Zalutskiy approach for the angular sector.

# Arguments
- `params`: QNMKerrCFParams object containing the Kerr parameters and initial guess
- `alg`: Nonlinear algorithm to use for the eigenvalue search (default: RobustMultiNewton with AutoFiniteDiff)
- `kwargs`: Additional keyword arguments to pass to the nonlinear solver
"""
function qnm_kerr_cf(
    params::QNMKerrCFParams{TR},
    alg::Union{AbstractNonlinearAlgorithm,Nothing}=RobustMultiNewton(;
        autodiff=AutoFiniteDiff()
    );
    kwargs...,
) where {TR<:AbstractFloat}
    @unpack_QNMKerrCFParams params

    cache = QNMKerrCFCache{TR}(params)
    ω0 = SA[real(ω_guess), imag(ω_guess)]
    prob = NonlinearProblem(qnm_kerr_cf_δ!, ω0, cache)
    sol = solve(prob, alg; kwargs...)
    ω = Complex{TR}(sol.u[1], sol.u[2])

    return ω
end

export qnm_kerr_radial_cf,
    QNMKerrCFParams, QNMKerrCFCache, qnm_kerr_cf_δ!, qnm_kerr_cf, @unpack_QNMKerrCFParams
