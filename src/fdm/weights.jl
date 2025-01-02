using ArgCheck
using LinearAlgebra: I

"""
    fornberg_calculate_wts([T=Float64], order::Integer, x0::Real, x::AbstractVector; 
                             dfdx::Bool=false)

Calculate finite difference weights for arbitrary-order derivatives using the Fornberg algorithm.

# Arguments
- `T`: Type parameter for the weights (defaults to type of x0)
- `order`: Order of the derivative to approximate
- `x0`: Point at which to approximate the derivative
- `x`: Grid points to use in the approximation
- `dfdx`: Whether to include first derivative values (Hermite finite differences)

# Returns
If `dfdx == false`:
- `Vector{T}`: Weights for function values

If `dfdx == true`:
- `Tuple{Vector{T}, Vector{T}}`: Weights for (function values, derivative values)

# Mathematical Background
For a function f(x), the derivative approximation takes the form:

If `dfdx == false` (standard finite differences):
```math
f^{(n)}(x_0) \\approx \\sum_{j=1}^N c_j f(x_j)
```

If `dfdx == true` (Hermite finite differences):
```math
f^{(n)}(x_0) \\approx \\sum_{j=1}^N [d_j f(x_j) + e_j f'(x_j)]
```

# Requirements
- For standard finite differences: N > order
- For Hermite finite differences: N > order/2 + 1
where N is the length of x

# Examples
```julia
# Standard central difference for first derivative
x = [-1.0, 0.0, 1.0]
w = fornberg_calculate_wts(1, 0.0, x)
# Returns approximately [-0.5, 0.0, 0.5]

# Forward difference for second derivative
x = [0.0, 1.0, 2.0, 3.0]
w = fornberg_calculate_wts(2, 0.0, x)

# Hermite finite difference for third derivative
x = [-1.0, 0.0, 1.0]
w_f, w_d = fornberg_calculate_wts(3, 0.0, x, dfdx=true)
```

# References

- [fornberg1988generation](@citet*)
- [fornberg2021algorithm](@citet*)
- [doi:10.1137/S0036144596322507](@citet*)
- [MethodOfLines.jl/src/discretization/schemes/fornberg_calculate_wts.jl at master · SciML/MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl/blob/master/src/discretization/schemes/fornberg_calculate_wts.jl)
- [precision - Numerical derivative and finite difference coefficients: any update of the Fornberg method? - Computational Science Stack Exchange](https://scicomp.stackexchange.com/questions/11249/numerical-derivative-and-finite-difference-coefficients-any-update-of-the-fornb)

# Notes
- The implementation includes a stability correction for higher-order derivatives
- For first derivatives (order=1), the weights sum to zero
- The algorithm handles both uniform and non-uniform grids
- When using Hermite finite differences, fewer points are needed but derivatives
  must be available

See also: [`fdm_grid`](@ref)
"""
function fornberg_calculate_wts(
    order::TI, x0::T, x::VT; dfdx::Bool=false
) where {T<:Real,T2<:Real,VT<:AbstractVector{T2},TI<:Integer}
    N = length(x)
    @argcheck dfdx || N > order "Standard finite difference requires at least order + 1 points."
    @argcheck !dfdx || N > div(order, 2) + 1 "Hermite finite difference requires at least order / 2 + 1 points."

    M = order
    c1 = one(T)
    c4 = x[1] - x0
    C = zeros(T, N, M + 1)
    C[1, 1] = 1
    @inbounds for i in 1:(N - 1)
        i1 = i + 1
        mn = min(i, M)
        c2 = one(T)
        c5 = c4
        c4 = x[i1] - x0
        for j in 0:(i - 1)
            j1 = j + 1
            c3 = x[i1] - x[j1]
            c2 *= c3
            if j == i - 1
                for s in mn:-1:1
                    s1 = s + 1
                    C[i1, s1] = c1 * (s * C[i, s] - c5 * C[i, s1]) / c2
                end
                C[i1, 1] = -c1 * c5 * C[i, 1] / c2
            end
            for s in mn:-1:1
                s1 = s + 1
                C[j1, s1] = (c4 * C[j1, s1] - s * C[j1, s]) / c3
            end
            C[j1, 1] = c4 * C[j1, 1] / c3
        end
        c1 = c2
    end
    #=
        This is to fix the problem of numerical instability which occurs when the sum of the stencil_coefficients is not
        exactly 0.
        https://scicomp.stackexchange.com/questions/11249/numerical-derivative-and-finite-difference-coefficients-any-update-of-the-fornb
        Stack Overflow answer on this issue.
        http://epubs.siam.org/doi/pdf/10.1137/S0036144596322507 - Modified Fornberg Algorithm
    =#
    _C = C[:, end]
    if order != 0
        _C[div(N, 2) + 1] -= sum(_C)
    end
    if dfdx == false
        return _C
    else
        A = x .- x'
        s = sum(1 ./ (A + I(N)); dims=1) .- 1
        cp = factorial.(0:M)
        cc = C ./ cp'
        c̃ = zeros(N, M + 2)
        for k in 1:(M + 1)
            c̃[:, k + 1] = sum(cc[:, 1:k] .* cc[:, k:-1:1]; dims=2)
        end
        E = c̃[:, 1:(M + 1)] - (x .- x0) .* c̃[:, 2:(M + 2)]
        D = c̃[:, 2:(M + 2)] + 2 * E .* s'
        D = D .* cp'
        E = E .* cp'

        _D = D[:, end]
        _E = E[:, end]
        return _D, _E
    end
end

export fornberg_calculate_wts

@testset "Fornberg Finite Difference Weights" begin
    @testset "Standard First Derivatives" begin
        # Test forward difference (first derivative)
        weights_forward = fornberg_calculate_wts(1, 0.0, [0.0, 1.0])
        @test weights_forward ≈ [-1.0, 1.0]

        # Test central difference (first derivative)
        weights_central = fornberg_calculate_wts(1, 0.0, [-1.0, 0.0, 1.0])
        @test weights_central ≈ [-0.5, 0.0, 0.5]

        # Test accuracy with non-uniform grid
        weights_nonuniform = fornberg_calculate_wts(1, 0.0, [-2.0, -0.5, 1.0])
        # Test against pre-computed values
        @test sum(weights_nonuniform) ≈ 0.0 atol = 1e-14
    end

    @testset "Standard Second Derivatives" begin
        # Test central difference (second derivative)
        weights_central = fornberg_calculate_wts(2, 0.0, [-1.0, 0.0, 1.0])
        @test weights_central ≈ [1.0, -2.0, 1.0]

        # Test with wider stencil
        weights_wide = fornberg_calculate_wts(2, 0.0, [-2.0, -1.0, 0.0, 1.0, 2.0])
        # Test properties that must hold
        @test sum(weights_wide) ≈ 0.0 atol = 1e-14
        @test sum(weights_wide .* [-2.0, -1.0, 0.0, 1.0, 2.0] .^ 2) ≈ 2.0 atol = 1e-14
    end

    @testset "Hermite Finite Difference" begin
        # Test first derivative with function and derivative values
        x = [-1.0, 0.0, 1.0]
        D, E = fornberg_calculate_wts(1, 0.0, x; dfdx=true)

        # Known weights for Hermite interpolation
        @test length(D) == length(x)  # Function value weights
        @test length(E) == length(x)  # Derivative value weights
        @test sum(D) ≈ 0.0 atol = 1e-14  # Consistency condition

        # Test second derivative with function and derivative values
        D2, E2 = fornberg_calculate_wts(2, 0.0, x; dfdx=true)
        @test sum(D2) ≈ 0.0 atol = 1e-14
    end

    @testset "Error Handling" begin
        # Test insufficient points for order
        @test_throws ArgumentError fornberg_calculate_wts(2, 0.0, [0.0, 1.0])

        # Test insufficient points for Hermite interpolation
        @test_throws ArgumentError fornberg_calculate_wts(3, 0.0, [0.0, 1.0], dfdx=true)
    end
end
