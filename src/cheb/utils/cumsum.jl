"""
    ChebCumsumOp{TR<:AbstractFloat,TI<:Integer}

A pre-allocated operator for computing the indefinite integral (cumulative sum) of a function
represented in the Chebyshev basis.

# Fields
- `n::TI`: Number of coefficients in the expansion
- `tmp::Vector{TR}`: Temporary storage for padded coefficients
- `result::Vector{TR}`: Storage for the result coefficients
- `v::Vector{TI}`: Pre-computed alternating signs [1, -1, 1, -1, ...]

# Type Parameters
- `TR`: The floating-point type for coefficients (e.g., Float64)
- `TI`: The integer type for indexing and signs (e.g., Int64)
"""
struct ChebCumsumOp{TR<:AbstractFloat,TI<:Integer}
    n::TI              # Number of coefficients
    tmp::Vector{TR}    # Temporary storage for padded coefficients
    result::Vector{TR} # Result storage
    v::Vector{TI}      # Pre-computed alternating signs

    function ChebCumsumOp(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
        # Pre-allocate workspace
        tmp = Vector{TR}(undef, n + 2)
        result = Vector{TR}(undef, n + 1)

        # Pre-compute alternating signs [1, -1, 1, -1, ...]
        v = ones(TI, n)
        @inbounds for i in 2:2:n
            v[i] = -one(TI)
        end

        return new{TR,TI}(n, tmp, result, v)
    end
end

"""
    (op::ChebCumsumOp)(f::AbstractVector)

Compute the indefinite integral of a function given its Chebyshev coefficients.

# Arguments
- `f`: Vector of Chebyshev coefficients of the function to be integrated

# Returns
- Vector of Chebyshev coefficients of the indefinite integral

# Notes
The integration constant is chosen such that f(-1) = 0. The computation follows
the recurrence relation for Chebyshev integration:
1. b₂ = c₁ - c₃/2
2. bᵣ = (cᵣ₋₁ - cᵣ₊₁)/(2r) for r > 1
3. b₀ is computed to ensure f(-1) = 0

where bᵢ are the coefficients of the integral and cᵢ are the input coefficients.
"""
function (op::ChebCumsumOp{TR,TI})(
    f::VT
) where {TR<:AbstractFloat,TI<:Integer,VT<:AbstractVector{TR}}
    @argcheck length(f) == op.n "length(f) must be equal to n"

    n = length(f)
    tmp = op.tmp
    result = op.result
    v = op.v

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

"""
    cheb_cumsum(f::AbstractVector)

Compute the indefinite integral of a function given its Chebyshev coefficients.
This is a convenience wrapper that creates a temporary ChebCumsumOp.

# Arguments
- `f`: Vector of Chebyshev coefficients of the function to be integrated

# Returns
- Vector of Chebyshev coefficients of the indefinite integral

# Example
```julia
# Integrate cos(x)
n = 15
f = cos.(cheb1_pts(n))
f_coeffs = cheb1_vals2coeffs(f)
If_coeffs = cheb_cumsum(f_coeffs)  # Coefficients of sin(x) - sin(-1)
```
"""
function cheb_cumsum(f::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    n = length(f)
    op = ChebCumsumOp(TR, n)
    return op(f)
end

export cheb_cumsum, ChebCumsumOp

@testset "cheb_cumsum" begin
    tol = 100 * eps()

    @testset "cheb1" begin
        n = 15
        f = cos.(cheb1_pts(n))
        f_coeffs = cheb1_vals2coeffs(f)
        If_coeffs = cheb_cumsum(f_coeffs)
        If = sin.(cheb1_pts(n)) .- sin(-1) # sin(x) - sin(-1) is the antiderivative of cos(x)
        If_coeffs_true = cheb1_vals2coeffs(If)
        @test norm(If_coeffs[1:(end - 1)] .- If_coeffs_true, Inf) < tol
    end

    @testset "cheb2" begin
        n = 15
        f = cos.(cheb2_pts(n))
        f_coeffs = cheb2_vals2coeffs(f)
        If_coeffs = cheb_cumsum(f_coeffs)
        If = sin.(cheb2_pts(n)) .- sin(-1) # sin(x) - sin(-1) is the antiderivative of cos(x)
        If_coeffs_true = cheb2_vals2coeffs(If)
        @test norm(If_coeffs[1:(end - 1)] .- If_coeffs_true, Inf) < tol
    end
end
