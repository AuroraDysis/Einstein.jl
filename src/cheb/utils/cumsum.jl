"""
    cumsum(f::Vector{TR}) where {TR<:AbstractFloat}

Compute the indefinite integral of a function represented in Chebyshev basis.

# Mathematical Background
Given a Chebyshev series of length n:
```math
G(x) = \\sum_{r=0}^{n-1} c_r T_r(x)
```
its integral is represented with a series of length n+1:
```math
\\int G(x)\\,dx = \\sum_{r=0}^{n} b_r T_r(x)
```

# Arguments
- `f::Vector{TR}`: Coefficients ``c_r`` of the Chebyshev series

# Returns
- Vector of coefficients ``b_r`` for the integral, with length n+1

# Mathematical Details
The coefficients are computed as follows:
```math
\\begin{align*}
b_0 &= \\sum_{r=1}^{n} (-1)^{r+1} b_r \\quad \\text{(constant of integration)} \\\\
b_1 &= c_0 - \\frac{c_2}{2} \\\\
b_r &= \\frac{c_{r-1} - c_{r+1}}{2r} \\quad \\text{for } r > 1
\\end{align*}
```
where ``c_{n+1} = c_{n+2} = 0``

The constant term ``b_0`` is chosen to make ``f(-1) = 0``.

# Examples
```julia
using Test

# Suppose we have 3 Chebyshev coefficients:
f = [1.0, 2.0, 3.0]

# Perform the continuous sum:
f_new = cumsum(f)

@test length(f_new) == length(f) + 1
println("Original coefficients: ", f)
println("New coefficients: ", f_new)
```

# References
- [mason2002chebyshev](@citet*), pp. 32-33.
- [chebfun/@chebtech/cumsum.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech/cumsum.m)

See also: [`diff`](@ref), [`sum`](@ref)
"""

struct ChebCumsumCache{TR<:AbstractFloat}
    tmp::Vector{TR}    # Temporary storage for padded coefficients
    result::Vector{TR} # Result storage
    v::Vector{TR}      # Pre-computed alternating signs

    function ChebCumsumCache{TR}(n::TI) where {TR<:AbstractFloat,TI<:Integer}
        # Pre-allocate workspace
        tmp = Vector{TR}(undef, n + 2)
        result = Vector{TR}(undef, n + 1)

        # Pre-compute alternating signs [1, -1, 1, -1, ...]
        v = ones(TI, n)
        @inbounds for i in 2:2:n
            v[i] = -one(TI)
        end

        return new{TR}(tmp, result, v)
    end
end

function cumsum(f::Vector{TR}) where {TR<:AbstractFloat}
    n = length(f)
    cache = ChebCumsumCache{TR}(n)
    return cumsum!(f, cache)
end

function cumsum!(f::Vector{TR}, cache::ChebCumsumCache{TR}) where {TR<:AbstractFloat}
    n = length(f)
    tmp = cache.tmp
    result = cache.result
    v = cache.v

    # Copy and pad input coefficients
    @inbounds begin
        tmp[1:n] .= f
        tmp[n + 1] = 0
        tmp[n + 2] = 0
    end

    # Compute interior coefficients
    @inbounds begin
        # b₂ = c₁ - c₃/2
        result[2] = tmp[1] - tmp[3] / 2

        # bᵣ = (cᵣ₋₁ - cᵣ₊₁)/(2r) for r > 1
        for r in 2:n
            result[r + 1] = (tmp[r] - tmp[r + 2]) / (2 * r)
        end
    end

    # Compute b₀ to ensure f(-1) = 0
    @inbounds begin
        result[1] = 0
        for i in 1:n
            result[1] += v[i] * result[i + 1]
        end
    end

    return result
end

export cumsum, cumsum!, ChebCumsumCache

@testset "cumsum" begin
    tol = 100 * eps()

    @testset "cheb1" begin
        n = 15
        f = cos.(cheb1_pts(n))
        f_coeffs = cheb1_vals2coeffs(f)
        If_coeffs = cumsum(f_coeffs)
        If = sin.(cheb1_pts(n)) .- sin(-1) # sin(x) - sin(-1) is the antiderivative of cos(x)
        If_coeffs_true = cheb1_vals2coeffs(If)
        @test norm(If_coeffs[1:(end - 1)] .- If_coeffs_true, Inf) < tol
    end

    @testset "cheb2" begin
        n = 15
        f = cos.(cheb2_pts(n))
        f_coeffs = cheb2_vals2coeffs(f)
        If_coeffs = cumsum(f_coeffs)
        If = sin.(cheb2_pts(n)) .- sin(-1) # sin(x) - sin(-1) is the antiderivative of cos(x)
        If_coeffs_true = cheb2_vals2coeffs(If)
        @test norm(If_coeffs[1:(end - 1)] .- If_coeffs_true, Inf) < tol
    end
end
