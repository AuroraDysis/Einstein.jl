"""
    cheb2_barywts([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}

Compute the barycentric weights for Chebyshev points of the 2nd kind.

# Arguments
- `TR`: Type parameter for the weights (e.g., Float64)
- `n`: Number of points

# References

- [chebfun/@chebtech2/barywts.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/%40chebtech2/barywts.m)

See also: [`bary`](@ref), [`cheb2_pts`](@ref)
"""
function cheb2_barywts(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    if n == 0
        return TR[]
    elseif n == 1
        return [one(TR)]
    end

    w = ones(TR, n)

    @inbounds begin
        half = one(TR) / 2
        w[(end - 1):-2:1] .= -1
        w[1] *= half
        w[end] = half
    end

    return w
end

function cheb2_barywts(n::TI) where {TI<:Integer}
    return cheb2_barywts(Float64, n)
end

export cheb2_barywts
