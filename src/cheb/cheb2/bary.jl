"""
    cheb2_bary_wts([TR=Float64], n::Integer)

Compute the barycentric weights for Chebyshev points of the second kind.

# Arguments
- `TR`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# Returns
- Vector of n barycentric weights

# Mathematical Details
For Chebyshev points of the second kind, the barycentric weights are:
```math
w_j = (-1)^j \\delta_j, \\quad j = 0,\\ldots,n-1
```
where ``\\delta_j`` is defined as:
```math
\\delta_j = \\begin{cases}
1/2 & j = 0 \\text{ or } j = n-1 \\\\
1 & \\text{otherwise}
\\end{cases}
```

These weights are optimized for numerical stability and efficiency in the barycentric
interpolation formula.

See also: [`bary`](@ref), [`cheb2_grid`](@ref)
"""
function cheb2_bary_wts(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    if n == 0
        return TR[]
    elseif n == 1
        return [one(TR)]
    end

    w = ones(TR, n)

    @inbounds begin
        w[(end - 1):-2:1] .= -1
        w[1] /= 2
        w[end] = one(TR) / 2
    end

    return w
end

function cheb2_bary_wts(n::TI) where {TI<:Integer}
    return cheb2_bary_wts(Float64, n)
end

export cheb2_bary_wts

@testset "cheb2 Barycentric Interpolation" begin
    @testset "Barycentric weights" begin
        # Test n=1 case
        w1 = cheb2_bary_wts(Float64, 1)
        @test w1 ≈ [1.0]

        # Test n=2 case
        w2 = cheb2_bary_wts(Float64, 2)
        @test w2 ≈ [-0.5, 0.5]

        # Test n=5 case
        w5 = cheb2_bary_wts(Float64, 5)
        @test w5 ≈ [0.5, -1.0, 1.0, -1.0, 0.5]
    end

    @testset "Polynomial interpolation" begin
        n = 5
        x = cheb2_grid(Float64, n)
        w = cheb2_bary_wts(Float64, n)

        # Test exact interpolation of quadratic polynomial
        f = @. 2x^2 - x + 1
        x0 = 0.3  # Test point
        exact = 2x0^2 - x0 + 1
        interp = bary(w, x, f, x0)
        @test interp ≈ exact rtol = 1e-14

        # Test interpolation at nodes (should be exact)
        for i in 1:n
            @test bary(w, x, f, x[i]) ≈ f[i] rtol = 1e-14
        end
    end

    @testset "Trigonometric function interpolation" begin
        n = 32  # More points for better accuracy
        x = cheb2_grid(Float64, n)
        w = cheb2_bary_wts(Float64, n)

        # Test interpolation of sin(πx)
        f = @. sin(π * x)
        test_points = [-0.7, -0.3, 0.0, 0.4, 0.8]
        f_exact = sin.(π * test_points)
        f_interp = [bary(w, x, f, x0) for x0 in test_points]
        @test f_interp ≈ f_exact rtol = 1e-12
    end
end
