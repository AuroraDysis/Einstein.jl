"""
MIT License

Copyright (c) 2025 Zhen Zhong
Copyright (c) 2022 SciML Open Source Scientific Machine Learning Organization

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

@doc raw"""
    fdm_weights_fornberg([TR=Float64], order::Integer, x0::Real, x::AbstractVector; 
                             hermite::Bool=false)

Calculate finite difference weights for arbitrary-order derivatives using the Fornberg algorithm.
Taken from [SciML/MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl).

# Arguments
- `TR`: Type parameter for the weights (defaults to type of x0)
- `order`: Order of the derivative to approximate
- `x0`: Point at which to approximate the derivative
- `x`: Grid points to use in the approximation
- `hermite`: Whether to include first derivative values (Hermite finite differences)

# Returns
If `hermite == false`:
- `Vector{TR}`: Weights for standard finite differences

If `hermite == true`:
- `Tuple{Vector{TR}, Vector{TR}}`: Weights for Hermite finite differences

# Requirements
- For standard finite differences: N > order
- For Hermite finite differences: N > order/2 + 1
where N is the length of x

# Examples
```julia
# Standard central difference for first derivative
x = [-1.0, 0.0, 1.0]
w = fdm_weights_fornberg(1, 0.0, x)
# Returns approximately [-0.5, 0.0, 0.5]

# Forward difference for second derivative
x = [0.0, 1.0, 2.0, 3.0]
w = fdm_weights_fornberg(2, 0.0, x)

# Hermite finite difference for third derivative
x = [-1.0, 0.0, 1.0]
w_f, w_d = fdm_weights_fornberg(3, 0.0, x, hermite=true)
```

# References

- [MethodOfLines.jl/src/discretization/schemes/fornberg_calculate_weights.jl at master · SciML/MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl/blob/master/src/discretization/schemes/fornberg_calculate_weights.jl)
- [fornberg1988generation](@citet*)
- [fornberg2021algorithm](@citet*)
- [Fornberg1998](@citet*)
- [precision - Numerical derivative and finite difference coefficients: any update of the Fornberg method? - Computational Science Stack Exchange](https://scicomp.stackexchange.com/questions/11249/numerical-derivative-and-finite-difference-coefficients-any-update-of-the-fornb)
"""
function fdm_weights_fornberg(
    order::Integer, x0::TR, x::AbstractVector{TR}; hermite::Bool=false
) where {TR<:Real}
    N = length(x)
    @argcheck hermite || N > order "Standard finite difference requires at least order + 1 points."
    @argcheck !hermite || N > div(order, 2) + 1 "Hermite finite difference requires at least order / 2 + 1 points."

    M = order
    c1 = one(TR)
    c4 = x[1] - x0
    C = zeros(TR, N, M + 1)
    C[1, 1] = 1
    @inbounds for i in 1:(N - 1)
        i1 = i + 1
        mn = min(i, M)
        c2 = one(TR)
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
    if hermite == false
        return _C
    else
        A = x .- x'
        s = sum(1 ./ (A + I(N)); dims=1) .- 1
        cp = factorial.(0:M)
        cc = C ./ cp'
        c̃ = zeros(TR, N, M + 2)
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

"""
    fdm_central_weights([TR=Rational{TI}], der_order::TI, acc_order::TI) where {TR<:Real, TI<:Integer}

Generate central finite difference coefficients for a given derivative and accuracy order.

# Arguments
- `der_order::Integer`: The order of the derivative to approximate
- `acc_order::Integer`: The desired order of accuracy (must be even)

# Returns
Vector of rational coefficients for the finite difference stencil
"""
function fdm_central_weights(
    ::Type{TR}, der_order::Integer, acc_order::Integer
) where {TR<:Real}
    @argcheck acc_order % 2 == 0 "Only even orders are supported for central FDM stencils."

    num_coeffs = fdm_central_size(der_order, acc_order)
    num_side = div(num_coeffs - 1, 2)
    local_grid = collect(TR, (-num_side):num_side)
    return fdm_weights_fornberg(der_order, zero(TR), local_grid)
end

function fdm_central_weights(der_order::TI, acc_order::TI) where {TI<:Integer}
    return fdm_central_weights(Rational{TI}, der_order, acc_order)
end

"""
    fdm_hermite_weights([TR=Rational{TI}], der_order::TI, acc_order::TI) where {TR<:Real, TI<:Integer}

Generate Hermite-type finite difference coefficients that include function value and derivative information.

# Arguments
- `der_order::Integer`: The order of the derivative to approximate (must be ≥ 2)
- `acc_order::Integer`: The desired order of accuracy
    * For der_order 2,3,6,7,10,11...: acc_order must be 4,8,12...
    * For der_order 4,5,8,9,12...: acc_order must be 2,6,10...

# Returns
Vector of rational coefficients for the Hermite-type finite difference stencil
"""
function fdm_hermite_weights(
    ::Type{TR}, der_order::Integer, acc_order::Integer
) where {TR<:Real}
    @argcheck der_order >= 2 "Only derivative order greater than or equal to 2 are supported for Hermite-type finite difference."

    if mod(div(der_order, 2), 2) == 1
        # acc_order must be 4,8,12... for der order 2,3,6,7,10,11...
        @argcheck acc_order % 4 == 0 "Only acc_order % 4 == 0 are supported for Hermite-type finite difference with der order 2,3,6,7,10,11..."
    else
        # acc_order must be 2,6,10... for der order 4,5,8,9,12...
        @argcheck acc_order % 4 == 2 "Only acc_order % 4 == 2 are supported for Hermite-type finite difference with der order 4,5,8,9,12..."
    end

    num_coeffs = fdm_hermite_size(der_order, acc_order)
    num_side = div(num_coeffs - 1, 2)
    local_grid = collect(TR, (-num_side):num_side)
    return fdm_weights_fornberg(der_order, zero(TR), local_grid; hermite=true)
end

function fdm_hermite_weights(der_order::TI, acc_order::TI) where {TI<:Integer}
    return fdm_hermite_weights(Rational{TI}, der_order, acc_order)
end

"""
    fdm_extrapwts_left(extrap_order::Int)

Generate weights for left-sided extrapolation of order `extrap_order`.

# Arguments
- `extrap_order::Int`: Order of extrapolation

# Returns
Vector of rational coefficients for left-sided extrapolation
"""
function fdm_extrapwts_left(extrap_order::Int)
    return fdm_weights_fornberg(0, 0, 1:extrap_order)
end

"""
    fdm_extrapwts_right(extrap_order::Int)

Generate weights for right-sided extrapolation of order `extrap_order`.

# Arguments
- `extrap_order::Int`: Order of extrapolation

# Returns
Vector of rational coefficients for right-sided extrapolation
"""
function fdm_extrapwts_right(extrap_order::Int)
    return fdm_weights_fornberg(0, 0, extrap_order:-1:1)
end

"""
    fdm_boundary_weights([TR=Rational{TI}], der_order::TI, acc_order::TI) where {TR<:Real, TI<:Integer}

Generate finite difference coefficients for shifted boundary conditions.

# Arguments
- `der_order::Integer`: The order of the derivative to approximate
- `acc_order::Integer`: The desired order of accuracy

# Returns
Tuple of left and right shifted boundary finite difference coefficients
The coefficients are stored in a matrix with the columns representing the different grid points.
The columns are ordered from the leftmost grid point to the rightmost grid point.
"""
function fdm_boundary_weights(
    ::Type{TR}, der_order::Integer, acc_order::Integer
) where {TR<:Real}
    num_coeffs = fdm_boundary_size(der_order, acc_order)
    num_central = fdm_central_size(der_order, acc_order)
    num_side = div(num_central - 1, 2)

    D_left = zeros(TR, num_coeffs, num_side)
    D_right = zeros(TR, num_coeffs, num_side)

    local_grid = zeros(TR, num_coeffs)
    @inbounds for i in 1:num_side
        for j in 1:num_coeffs
            local_grid[j] = -(i - 1) + (j - 1)
        end
        D_left[:, i] = fdm_weights_fornberg(der_order, zero(TR), local_grid)

        for j in 1:num_coeffs
            local_grid[end - j + 1] = (i - 1) - (j - 1)
        end
        D_right[:, end - i + 1] = fdm_weights_fornberg(der_order, zero(TR), local_grid)
    end

    return D_left, D_right
end

function fdm_boundary_weights(der_order::TI, acc_order::TI) where {TI<:Integer}
    return fdm_boundary_weights(Rational{TI}, der_order, acc_order)
end

export fdm_weights_fornberg,
    fdm_central_weights,
    fdm_hermite_weights,
    fdm_extrapwts_right,
    fdm_extrapwts_left,
    fdm_boundary_weights
