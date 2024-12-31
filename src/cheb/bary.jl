"""
    cheb2_interp_wts(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}

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

See also: [`cheb2_interp`](@ref), [`cheb2_grid`](@ref)
"""
function cheb2_interp_wts(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    w = ones(TR, n)

    if n == 1
        return w
    end

    @inbounds begin
        w[(end - 1):-2:1] .= -1
        w[1] /= 2
        w[end] = one(TR) / 2
    end

    return w
end

"""
    cheb2_interp(w::VT1, x::VT2, f::VT3, x0::TR) where {
        TR<:AbstractFloat,
        VT1<:AbstractVector{TR},
        VT2<:AbstractVector{TR},
        VT3<:AbstractVector{TR},
    }

Evaluate a polynomial interpolant using the barycentric interpolation formula [salzer1972lagrangian](@cite).

# Arguments
- `w`: Vector of barycentric weights
- `x`: Vector of interpolation points (typically Chebyshev points)
- `f`: Vector of function values at interpolation points
- `x0`: Point at which to evaluate the interpolant

# Returns
- Interpolated value at x0

# Mathematical Details
The barycentric interpolation formula is:
```math
p(x) = \\begin{cases}
f_j & \\text{if } x = x_j \\text{ for some } j \\\\
\\frac{\\sum_{j=0}^{n-1} \\frac{w_j}{x-x_j}f_j}{\\sum_{j=0}^{n-1} \\frac{w_j}{x-x_j}} & \\text{otherwise}
\\end{cases}
```

This formula provides a numerically stable way to evaluate the Lagrange interpolation
polynomial. When used with Chebyshev points and their corresponding barycentric weights,
it gives optimal interpolation properties.

# Examples
```julia
# Set up interpolation points and weights
n = 10
x = cheb2_grid(Float64, n)
w = cheb2_interp_wts(Float64, n)

# Function to interpolate
f = sin.(π .* x)

# Evaluate interpolant at a point
x0 = 0.5
y = cheb2_interp(w, x, f, x0)
```

# `chebfun` Reference

[chebfun/@chebtech2/bary.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/bary.m)

See also: [`cheb2_interp_wts`](@ref), [`cheb2_grid`](@ref)
"""
function cheb2_interp(
    w::VT1, x::VT2, f::VT3, x0::TR
) where {
    TR<:AbstractFloat,
    VT1<:AbstractVector{TR},
    VT2<:AbstractVector{TR},
    VT3<:AbstractVector{TR},
}
    p = zero(TR)
    q = zero(TR)

    @inbounds for i in eachindex(x)
        if x0 == x[i]
            return f[i]
        end

        wi = w[i] / (x0 - x[i])
        p += wi * f[i]
        q += wi
    end

    return p / q
end

export cheb2_interp_wts, cheb2_interp
