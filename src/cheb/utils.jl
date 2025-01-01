"""
    bary(w::VT1, x::VT2, f::VT3, x0::TR) where {
        TR<:AbstractFloat,
        VT1<:AbstractVector{TR},
        VT2<:AbstractVector{TR},
        VT3<:AbstractVector{TR},
    }

Evaluate a polynomial interpolant using the barycentric interpolation formula.

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
x = cheb2_pts(Float64, n)
w = cheb2_barywts(Float64, n)

# Function to interpolate
f = sin.(π .* x)

# Evaluate interpolant at a point
x0 = 0.5
y = bary(w, x, f, x0)
```

# Reference

- [salzer1972lagrangian](@citet*)
- [doi:10.1137/1.9781611975949](@citet*)
- [chebfun/@chebtech2/bary.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/bary.m)

See also: [`cheb1_barywts`](@ref), [`cheb1_pts`](@ref), [`cheb2_barywts`](@ref), [`cheb2_pts`](@ref)
"""
function bary(
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

export bary
