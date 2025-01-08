"""
    fdm_integrate_simpson(f::AbstractVector{T}, dx::T) where {T<:AbstractFloat}

Integrate a function `f` using Simpson's rule, given the grid spacing `dx`.
"""
function fdm_integrate_simpson(f::AbstractVector{T}, dx::T) where {T<:AbstractFloat}
    @inbounds retval =
        (
            17 * (f[1] + f[end]) +
            59 * (f[2] + f[end - 1]) +
            43 * (f[3] + f[end - 2]) +
            49 * (f[4] + f[end - 3])
        ) / 48
    @inbounds @simd for i in 5:(length(f) - 4)
        retval += f[i]
    end
    @inbounds return retval * dx
end

export fdm_integrate_simpson
