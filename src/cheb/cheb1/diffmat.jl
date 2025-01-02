"""
    cheb1_diffmat(n::TI, k::TI=1) where {TI<:Integer}

Construct a Chebyshev differentiation matrix for points of the first kind.

# Description
Creates a matrix that maps function values at `n` Chebyshev points of the first kind 
to values of the `k`-th derivative of the interpolating polynomial at those points.

# Arguments
- `n::Integer`: Number of Chebyshev points
- `k::Integer=1`: Order of the derivative (default: 1)

# Returns
- Matrix that computes the `k`-th derivative when multiplied with function values

# Examples
```julia
# First derivative matrix for 5 points
D1 = cheb1_diffmat(5)

# Second derivative matrix for 5 points
D2 = cheb1_diffmat(5, 2)

# Apply to function values
vals = [1.0, 2.0, 3.0, 4.0, 5.0]
deriv_vals = D1 * vals
```

# References
- Based on the MATLAB Chebfun implementation
"""
function cheb1_diffmat(::Type{TR}, n::TI, k::TI=1) where {TR<:AbstractFloat,TI<:Integer}
    x = cheb1_pts(TR, n)               # First kind points.
    w = cheb1_barywts(TR, n)           # Barycentric weights.
    t = cheb1_angles(TR, n)            # acos(x).
    D = bary_diffmat(x, w, k, t)       # Construct matrix.
    return D
end
