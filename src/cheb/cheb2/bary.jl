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
