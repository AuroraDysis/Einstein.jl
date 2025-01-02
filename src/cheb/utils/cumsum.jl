"""
    cheb_cumsum(f::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}
    ChebCumsumOp([TR=Float64], n::TI)(f::VT) where {TR<:AbstractFloat,TI<:Integer,VT<:AbstractVector{TR}}

Compute the indefinite integral of a function given its Chebyshev coefficients.

# Arguments
- `f`: Vector of Chebyshev coefficients of the function to be integrated

# References
- [chebfun/@chebtech/cumsum.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech/cumsum.m)
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

    function ChebCumsumOp(n::TI) where {TI<:Integer}
        return ChebCumsumOp{Float64,TI}(Float64, n)
    end
end

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
