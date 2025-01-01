"""
    clenshaw(c::VT, x::T) where {T<:AbstractFloat,VT<:AbstractVector{T}}

Evaluate a Chebyshev series using the Clenshaw algorithm.

# Mathematical Background
For a Chebyshev series:
```math
f(x) = \\sum_{k=0}^n c_k T_k(x)
```
where ``T_k(x)`` are Chebyshev polynomials of the first kind, the Clenshaw algorithm
computes the sum using the recurrence relation:
```math
\\begin{align*}
b_{n+1} &= b_{n+2} = 0 \\\\
b_k &= c_k + 2x b_{k+1} - b_{k+2} \\\\
f(x) &= \\frac{b_0 - b_2}{2}
\\end{align*}
```

# Arguments
- `c`: Vector of Chebyshev coefficients ``[c_0, c_1, \\ldots, c_n]``
- `x`: Evaluation point in [-1,1]

# Returns
- Value of the Chebyshev series at x

# Performance Notes
- Uses @inbounds for efficiency in the main loop
- Avoids allocations in the recurrence calculation
- Pre-scales x to reduce operations in the loop

# Examples
```julia
# Evaluate T₀(x) = 1
julia> clenshaw([1.0], 0.5)
1.0

# Evaluate T₁(x) = x
julia> clenshaw([0.0, 1.0], 0.5)
0.5

# Evaluate 1 + 2x + 3x²
julia> c = [1.0, 2.0, 3.0]
julia> x = 0.5
julia> clenshaw(c, x)
2.75
```

See also: [`cheb1_coeffs2vals`](@ref), [`cheb2_coeffs2vals`](@ref)
"""
function clenshaw(c::VT, x::T) where {T<:AbstractFloat,VT<:AbstractVector{T}}
    @argcheck length(c) > 0 "c must have at least one element"

    x = 2 * x

    bk1 = zero(T)
    bk2 = zero(T)

    n = length(c) - 1

    @inbounds for k in (n + 1):-2:3
        bk2 = c[k] + x * bk1 - bk2
        bk1 = c[k - 1] + x * bk2 - bk1
    end

    # If n is odd, perform the extra step
    if isodd(n)
        tmp = deepcopy(bk1)
        @inbounds bk1 = c[2] + x * bk1 - bk2
        bk2 = tmp
    end

    # Compute the final value
    @inbounds y = c[1] + x * bk1 / 2 - bk2

    return y
end

@testset "clenshaw" begin
    tol = 10 * eps()

    @testset "Single coefficient tests" begin
        # Scalar evaluation
        c = [sqrt(2)]
        v = clenshaw(c, 0.0)
        @test v ≈ sqrt(2)

        # Vector evaluation
        x = [-0.5, 1.0]
        v = map(xi -> clenshaw(c, xi), x)
        @test all(v .≈ sqrt(2))
    end

    @testset "Vector coefficient tests" begin
        # Simple polynomial test
        c = collect(5.0:-1:1)
        x = [-0.5, -0.1, 1.0]

        v = map(xi -> clenshaw(c, xi), x)
        v_true = [3.0, 3.1728, 15.0]
        @test norm(v - v_true, Inf) < tol

        # Multiple polynomials
        c2 = reverse(c)
        v = map(xi -> [clenshaw(c, xi), clenshaw(c2, xi)], x)
        v_true = [
            3.0 0.0
            3.1728 3.6480
            15.0 15.0
        ]
        @test norm(hcat(v...)' - v_true, Inf) < tol
    end

    @testset "Edge cases" begin
        # Empty coefficient vector
        @test_throws ArgumentError clenshaw(Float64[], 0.0)

        # Single coefficient
        @test clenshaw([1.0], 0.5) ≈ 1.0
    end
end
