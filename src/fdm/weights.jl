@doc raw"""
    fdm_fornbergwts([T=Float64], order::Integer, x0::Real, x::AbstractVector; 
                             dfdx::Bool=false)

Calculate finite difference weights for arbitrary-order derivatives using the Fornberg algorithm.
Taken from [SciML/MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl).

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
f^{(n)}(x_0) \approx \sum_{j=1}^N c_j f(x_j)
```

If `dfdx == true` (Hermite finite differences):
```math
f^{(n)}(x_0) \approx \sum_{j=1}^N [d_j f(x_j) + e_j f'(x_j)]
```

# Requirements
- For standard finite differences: N > order
- For Hermite finite differences: N > order/2 + 1
where N is the length of x

# Examples
```julia
# Standard central difference for first derivative
x = [-1.0, 0.0, 1.0]
w = fdm_fornbergwts(1, 0.0, x)
# Returns approximately [-0.5, 0.0, 0.5]

# Forward difference for second derivative
x = [0.0, 1.0, 2.0, 3.0]
w = fdm_fornbergwts(2, 0.0, x)

# Hermite finite difference for third derivative
x = [-1.0, 0.0, 1.0]
w_f, w_d = fdm_fornbergwts(3, 0.0, x, dfdx=true)
```

# References

- [MethodOfLines.jl/src/discretization/schemes/fdm_fornbergwts.jl at master · SciML/MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl/blob/master/src/discretization/schemes/fdm_fornbergwts.jl)
- [fornberg1988generation](@citet*)
- [fornberg2021algorithm](@citet*)
- [doi:10.1137/S0036144596322507](@citet*)
- [precision - Numerical derivative and finite difference coefficients: any update of the Fornberg method? - Computational Science Stack Exchange](https://scicomp.stackexchange.com/questions/11249/numerical-derivative-and-finite-difference-coefficients-any-update-of-the-fornb)
"""
function fdm_fornbergwts(
    order::Integer, x0::T, x::AbstractVector{T}; dfdx::Bool=false
) where {T<:Real}
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
        c̃ = zeros(T, N, M + 2)
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

export fdm_fornbergwts
