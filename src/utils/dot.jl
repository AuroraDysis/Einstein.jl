"""
    dot_xsum(x::StridedVector{T}, y::StridedVector{T}) where {T<:Real}

Compute the dot product of two vectors using extended precision accumulation.
Uses the `xsum` package for improved numerical accuracy.

# Note
- Very fast for large vectors, but a bit slower than Kahan summation for small vectors.
- Both input vectors must have the same length, which is not checked for performance reasons.

# References
- [neal2015fastexactsummationusing](@citet*)
- [JuliaMath/Xsum.jl](https://github.com/JuliaMath/Xsum.jl)
- [Radford Neal / xsum · GitLab](https://gitlab.com/radfordneal/xsum)
"""
function dot_xsum(x::StridedVector{Float64}, y::StridedVector{Float64})
    acc = XAccumulator(0.0)
    @inbounds for i in eachindex(x)
        accumulate!(acc, x[i] * y[i])
    end
    return float(acc)
end

"""
    dot_kahan(v1::StridedVector{T}, v2::StridedVector{T}) where {T<:Number}

Compute the dot product using Kahan summation algorithm to reduce numerical errors.

# Note
- Slower than `dot_xsum` for large vectors, but faster for small vectors.
- Similar performance to `dot_kahan_opt`
- Both input vectors must have the same length, which is not checked for performance reasons.
"""
function dot_kahan(v1::StridedVector{T}, v2::StridedVector{T}) where {T<:Number}
    s = zero(T)
    c = zero(T)
    y = zero(T)
    t = zero(T)
    j = 1
    n = length(v1)

    @inbounds for j in 1:n
        vj = v1[j] * v2[j]
        y = vj - c
        t = s
        s += y
        c = (s - t) - y
    end

    return s
end

"""
    dot_kahan_opt(v1::StridedVector{T}, v2::StridedVector{T}) where {T<:Number}

Optimized version of Kahan summation for dot product computation.

# Note
- Slower than `dot_xsum` for large vectors, but faster for small vectors.
- Similar performance to `dot_kahan`
- Uses loop unrolling for better performance while maintaining Kahan summation's
numerical stability. Processes two elements per iteration when possible.

# References
- [Radford Neal / xsum · GitLab](https://gitlab.com/radfordneal/xsum)
"""
function dot_kahan_opt(v1::StridedVector{T}, v2::StridedVector{T}) where {T<:Number}
    s = zero(T)
    c = zero(T)
    y = zero(T)
    t = zero(T)
    j = 2
    n = length(v1)

    @inbounds while j <= n
        vjm1 = v1[j - 1] * v2[j - 1]
        y = vjm1 - c
        t = s
        s += y
        c = (s - t) - y
        vj = v1[j] * v2[j]
        y = vj - c
        t = s
        s += y
        c = (s - t) - y
        j += 2
    end

    @inbounds for k in (j - 1):n
        vk = v1[k] * v2[k]
        y = vk - c
        t = s
        s += y
        c = (s - t) - y
    end

    return s
end

export dot_xsum, dot_kahan, dot_kahan_opt
