"""
    fdm_integrate_simpson(f::AbstractVector{T}, dx::T) where {T<:AbstractFloat}

Integrate a function `f` using Simpson's rule, given the grid spacing `dx`.
"""
function fdm_integrate_simpson(f::AbstractVector{T}, dx::T) where {T<:AbstractFloat}
    @argcheck length(f) >= 4 "f must have at least 4 elements"

    @inbounds val =
        (
            17 * (f[1] + f[end]) +
            59 * (f[2] + f[end - 1]) +
            43 * (f[3] + f[end - 2]) +
            49 * (f[4] + f[end - 3])
        ) / 48
    @inbounds @simd for i in 5:(length(f) - 4)
        val += f[i]
    end
    @inbounds return val * dx
end

function fdm_integrate_trapezoidal(f::AbstractVector{T}, dx::T) where {T<:AbstractFloat}
    @argcheck length(f) >= 3 "f must have at least 3 elements"

    @inbounds val = f[2]
    @inbounds @simd for i in 3:(length(f) - 1)
        val += f[i]
    end
    @inbounds return dx * (val + (f[1] + f[end]) / 2)
end

export fdm_integrate_simpson, fdm_integrate_trapezoidal
