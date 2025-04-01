"""
    sum_xsum(vec::AbstractVector{Float64})

Compute the sum of a vector using extended precision accumulation.
Uses the `xsum` package for improved numerical accuracy.

# Note
- Very fast for large vectors, but a bit slower than Kahan summation for small vectors (n ⪅ 80)

# References
- [neal2015fastexactsummationusing](@citet*)
- [JuliaMath/Xsum.jl](https://github.com/JuliaMath/Xsum.jl)
- [Radford Neal / xsum · GitLab](https://gitlab.com/radfordneal/xsum)
"""
@inline function sum_xsum(vec::AbstractVector{Float64})
    return xsum(vec)
end

"""
    sum_kahan(v::AbstractVector{T}) where {T<:Number}

Neumaier's variant of Kahan summation algorithm to reduce numerical errors.

# Note
- Slower than `sum_xsum` for large vectors, but faster for small vectors.
- Uses loop unrolling for better performance while maintaining Kahan summation's
numerical stability. Processes two elements per iteration when possible.

# References
- [JuliaMath/KahanSummation.jl](https://github.com/JuliaMath/KahanSummation.jl)
"""
function sum_kahan(v::AbstractVector{T}) where {T<:Number}
    @inbounds begin
        it = iterate(v)
        c = zero(T)

        # return 0 if empty(v)
        if it === nothing
            return c
        end

        vi, i = it
        s = vi - c

        while (it = iterate(v, i)) !== nothing
            vi, i = it
            t = s + vi
            if abs(s) >= abs(vi)
                c -= ((s - t) + vi)
            else
                c -= ((vi - t) + s)
            end
            s = t
        end

        return s - c
    end
end

export sum_xsum, sum_kahan
