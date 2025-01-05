@enumx SchwPotential ReggeWheeler Zerilli
@enumx BoundaryCondition Natural Dirichlet

"""
    qnm_schwpep(::Type{TR}, s::Integer, ℓ::Integer, cheb_n::Integer, potential::SchwPotential.T; σ_min::TR=zero(TR), σ_max::TR=one(TR), lo_bc::BoundaryCondition.T=Natural, hi_bc::BoundaryCondition.T=Natural)

Construct a polynomial eigenvalue problem for the Schwarzschild spacetime using the Ultraspherical spectral method.

# Arguments
- `s::Integer`: Spin
- `ℓ::Integer`: Angular number
- `cheb_n::Integer`: Number of Chebyshev points
- `potential::SchwPotential.T`: Potential type
- `σ_min::TR=zero(TR)`: Minimum value of the radial coordinate (hyperboloidal slicing)
- `σ_max::TR=one(TR)`: Maximum value of the radial coordinate (hyperboloidal slicing)
- `lo_bc::BoundaryCondition.T=Natural`: Boundary condition at the lower boundary, either `Natural` or `Dirichlet`.
- `hi_bc::BoundaryCondition.T=Natural`: Boundary condition at the upper boundary, either `Natural` or `Dirichlet`.

# Returns
Polynomial eigenvalue problem

# References
- [Jaramillo:2020tuu](@citet*)
"""
function qnm_schwpep(
    ::Type{TR},
    s::Integer,
    ℓ::Integer,
    cheb_n::Integer,
    potential::SchwPotential.T;
    σ_min::TR=zero(TR),
    σ_max::TR=one(TR),
    lo_bc::BoundaryCondition.T=Natural,
    hi_bc::BoundaryCondition.T=Natural,
) where {TR<:AbstractFloat,TI<:Integer}
    # Zerilli must have s = 2
    if potential == Zerilli
        @argcheck s == 2 "s must be 2 for Zerilli potential"
    end

    @argcheck ℓ >= s "ℓ must be greater than or equal to s"
    @argcheck cheb_n >= 2 "cheb_n must be greater than or equal to 2"

    chebSpace = Chebyshev(σ_min .. σ_max)
    ultraSpace = Ultraspherical(2, σ_min .. σ_max)
    conversion = Conversion(chebSpace, ultraSpace)
    conversionA1 = Conversion(Ultraspherical(1, σ_min .. σ_max), ultraSpace)

    σ = Fun(chebSpace)
    c02 = -(σ - 1) * σ^2
    c01 = (2 - 3 * σ) * σ
    c00 = if potential == ReggeWheeler
        -(ℓ * (ℓ + 1)) + (s^2 - 1) * σ
    else
        let n = TR(ℓ - 1) * TR(ℓ + 2) / 2
            -σ - 2n * (4n * σ + 4n * (n + 1) + 3σ^2) / (2n + 3σ)^2
        end
    end
    A0 = (c02 * 𝒟^2 + c01 * 𝒟 + c00):chebSpace

    c11 = im * (-1 + 2 * σ^2)
    c10 = 2im * σ
    A1 = (c11 * 𝒟 + c10):chebSpace
    A1c = conversionA1 * A1

    c20 = σ + 1
    A2 = (c20 * conversion):chebSpace

    A0m = zeros(TR, cheb_n, cheb_n)
    A1m = zeros(Complex{TR}, cheb_n, cheb_n)
    A2m = zeros(TR, cheb_n, cheb_n)
    A0m .= A0[1:cheb_n, 1:cheb_n]
    A1m .= A1c[1:cheb_n, 1:cheb_n]
    A2m .= A2[1:cheb_n, 1:cheb_n]

    if lo_bc == Dirichlet && hi_bc == Dirichlet
        # not test yet
        A0m[end - 1, 1:2:end] .= 1
        A0m[end - 1, 2:2:end] .= -1
        A0m[end, :] .= 1
        A1m[(end - 1):end, :] .= 0
        A2m[(end - 1):end, :] .= 0
    elseif lo_bc == Dirichlet && hi_bc == Natural
        # not test yet
        A0m[end, 1:2:end] .= 1
        A0m[end, 2:2:end] .= -1
        A1m[end, :] .= 0
        A2m[end, :] .= 0
    elseif lo_bc == Natural && hi_bc == Dirichlet
        A0m[end, :] .= 1
        A1m[end, :] .= 0
        A2m[end, :] .= 0
    end

    nep = PEP([A0m, A1m, A2m])

    return nep
end

export SchwPotential, BoundaryCondition, qnm_schwpep
