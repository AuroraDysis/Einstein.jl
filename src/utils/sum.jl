"""
    sum_xsum(vec::StridedVector{T}) where {T<:Union{Float32,Float64}}

Compute the sum of a vector using extended precision accumulation.
Uses the `xsum` package for improved numerical accuracy.

# Note
- Very fast for large vectors, but a bit slower than Kahan summation for small vectors (n ⪅ 80)

# References
- [neal2015fastexactsummationusing](@citet*)
- [JuliaMath/Xsum.jl](https://github.com/JuliaMath/Xsum.jl)
- [Radford Neal / xsum · GitLab](https://gitlab.com/radfordneal/xsum)
"""
@inline function sum_xsum(vec::StridedVector{Float64})
    return xsum(vec)
end

"""
    sum_kahan(v::StridedVector{T}) where {T<:Number}

Compute the sum using Kahan summation algorithm to reduce numerical errors.

# Note
- Slower than `sum_xsum` for large vectors, but faster for small vectors.
- Similar performance to `sum_kahan_opt`
"""
function sum_kahan(v::StridedVector{T}) where {T<:Number}
    s = zero(T)
    c = zero(T)
    y = zero(T)
    t = zero(T)
    j = 1
    n = length(v)

    @inbounds for j in 1:n
        y = v[j] - c
        t = s
        s += y
        c = (s - t) - y
    end

    return s
end

"""
    sum_kahan_opt(v::StridedVector{T}) where {T<:Number}

Optimized version of Kahan summation for vector sum computation.

# Note
- Slower than `sum_xsum` for large vectors, but faster for small vectors.
- Similar performance to `sum_kahan`
- Uses loop unrolling for better performance while maintaining Kahan summation's
numerical stability. Processes two elements per iteration when possible.

# References
- [Radford Neal / xsum · GitLab](https://gitlab.com/radfordneal/xsum)
"""
function sum_kahan_opt(v::StridedVector{T}) where {T<:Number}
    s = zero(T)
    c = zero(T)
    y = zero(T)
    t = zero(T)
    j = 2
    n = length(v)

    @inbounds while j <= n
        y = v[j - 1] - c
        t = s
        s += y
        c = (s - t) - y
        y = v[j] - c
        t = s
        s += y
        c = (s - t) - y
        j += 2
    end

    @inbounds for k in (j - 1):n
        y = v[k] - c
        t = s
        s += y
        c = (s - t) - y
    end

    return s
end

export sum_xsum, sum_kahan, sum_kahan_opt
