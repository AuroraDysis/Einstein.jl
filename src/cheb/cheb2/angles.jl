@doc raw"""
    cheb2_angles([TF=Float64], n::Integer) where {TF<:AbstractFloat}
    cheb2_angles!(θ::Vector{TF}, n::Integer) where {TF<:AbstractFloat}

Compute angles for Chebyshev points of the 2nd kind:
```math
\theta_k = \frac{k\pi}{n-1}, \quad k = n-1,\ldots,0
```

# Arguments
- `TF`: Type parameter for the angles (e.g., Float64)
- `n`: Number of points
"""
function cheb2_angles(::Type{TF}, n::Integer) where {TF<:AbstractFloat}
    @argcheck n >= 0 "n must be nonnegative"

    θ = Array{TF}(undef, n)
    cheb2_angles!(θ, n)

    return θ
end

function cheb2_angles(n::Integer)
    return cheb2_angles(Float64, n)
end

function cheb2_angles!(θ::Vector{TF}, n::Integer) where {TF<:AbstractFloat}
    @argcheck length(θ) == n "length(θ) must be equal to n"

    if n == 1
        θ[1] = convert(TF, π) / 2
        return nothing
    end

    nm1 = n - 1
    pi_over_nm1 = convert(TF, π) / nm1

    @inbounds for k in 0:nm1
        θ[n - k] = k * pi_over_nm1
    end

    return nothing
end

function angles(grid::ChebyshevGrid{TF,ChebyshevSecondKindNode}) where {TF}
    if !haskey(grid.cached, :angles)
        resize!(grid.angles, grid.n)
        cheb2_angles!(grid.angles, grid.n)
        grid.cached[:angles] = true
    end
    return grid.angles
end

export cheb2_angles, cheb2_angles!, angles
