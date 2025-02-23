@doc raw"""
    qnm_schw_cheb(::Type{TR}, s::Integer, ℓ::Integer, n::Integer, potential::SchwPType.T; σ_min::TR=zero(TR), σ_max::TR=one(TR), lo_bc::BCType.T=Natural, hi_bc::BCType.T=Natural)

Construct a polynomial eigenvalue problem for the Schwarzschild spacetime using the hyperboloidal coordinates and the ultraspherical spectral method.
The coordinate transformation from standard Schwarzschild coordinates to hyperboloidal coordinates is given by (we use $M = 1$ in the code):
```math
\begin{align}
\sigma &= \frac{2 M}{r} \\
\tau &= t + 2 M \left( \ln\sigma + \ln(1 - \sigma) - \frac{1}{\sigma} \right) \\
r_* &= 2 M \left(\frac{1}{\sigma} + \ln(1 - \sigma) - \ln\sigma \right) 
\end{align}
```

# Arguments
- `s::Integer`: Spin
- `ℓ::Integer`: Angular number
- `n::Integer`: Number of Chebyshev points
- `potential::SchwPType.T`: Potential type
- `σ_min::TR=zero(TR)`: Minimum value of the radial coordinate (hyperboloidal slicing)
- `σ_max::TR=one(TR)`: Maximum value of the radial coordinate (hyperboloidal slicing)
- `lo_bc::BCType.T=BCType.Natural`: Boundary condition at the lower boundary, either `Natural` or `Dirichlet`.
- `hi_bc::BCType.T=BCType.Natural`: Boundary condition at the upper boundary, either `Natural` or `Dirichlet`.

# Returns
A array of matrices representing the Polynomial eigenvalue problem, which can be solved using `qnm_polyeig` or solvers from the [NonlinearEigenproblems.jl](https://github.com/nep-pack/NonlinearEigenproblems.jl) package, such as `polyeig`.

# References
- [Jaramillo:2020tuu, PanossoMacedo:2023qzp](@citet*)
"""
function qnm_schw_cheb(
    ::Type{TR},
    s::Integer,
    ℓ::Integer,
    n::Integer,
    potential::SchwPType.T;
    σ_min::TR=zero(TR),
    σ_max::TR=one(TR),
    lo_bc::BCType.T=BCType.Natural,
    hi_bc::BCType.T=BCType.Natural,
) where {TR<:AbstractFloat}
    # Zerilli must have s = 2
    if potential == SchwPType.Zerilli
        @argcheck s == 2 "s must be 2 for Zerilli potential"
    end

    @argcheck ℓ >= s "ℓ must be greater than or equal to s"
    @argcheck n >= 2 "n must be greater than or equal to 2"

    dom = σ_min .. σ_max
    chebSpace = Chebyshev(dom)
    ultraSpace = Ultraspherical(2, dom)
    conversion = Conversion(chebSpace, ultraSpace)
    conversionA1 = Conversion(Ultraspherical(1, dom), ultraSpace)

    σ = Fun(chebSpace)
    c02 = -(σ - 1) * σ^2 / 16
    c01 = (2 - 3 * σ) * σ / 16
    c00 = if potential == SchwPType.ReggeWheeler
        (-(ℓ * (ℓ + 1)) + (s^2 - 1) * σ) / 16
    else
        let n = TR(ℓ - 1) * TR(ℓ + 2) / 2
            (-σ - 2n * (4n * σ + 4n * (n + 1) + 3σ^2) / (2n + 3σ)^2) / 16
        end
    end
    A0 = (c02 * 𝒟^2 + c01 * 𝒟 + c00):chebSpace

    c11 = 1im * (2 * σ^2 - 1) / 4
    c10 = 1im * σ / 2
    A1 = (c11 * 𝒟 + c10):chebSpace
    A1c = conversionA1 * A1

    c20 = σ + 1
    A2 = (c20 * conversion):chebSpace

    A0m = @view(A0[1:n, 1:n]) |> Matrix
    A1m = @view(A1c[1:n, 1:n]) |> Matrix
    A2m = @view(A2[1:n, 1:n]) |> Matrix

    if lo_bc == BCType.Dirichlet && hi_bc == BCType.Dirichlet
        # TODO: add tests for this case
        A0m[end - 1, 1:2:end] .= 1
        A0m[end - 1, 2:2:end] .= -1
        A0m[end, :] .= 1
        A1m[(end - 1):end, :] .= 0
        A2m[(end - 1):end, :] .= 0
    elseif lo_bc == BCType.Dirichlet && hi_bc == BCType.Natural
        # TODO: add tests for this case
        A0m[end, 1:2:end] .= 1
        A0m[end, 2:2:end] .= -1
        A1m[end, :] .= 0
        A2m[end, :] .= 0
    elseif lo_bc == BCType.Natural && hi_bc == BCType.Dirichlet
        # TODO: add tests for this case
        A0m[end, :] .= 1
        A1m[end, :] .= 0
        A2m[end, :] .= 0
    end

    nep = [A0m, A1m, A2m]

    return nep
end

export qnm_schw_cheb
