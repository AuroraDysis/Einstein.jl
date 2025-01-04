using Xsum

# TODO: benchmark this against the naive implementation
# The size of the vectors should be same
function dot_xsum(x::StridedVector{T}, y::StridedVector{T}) where {T<:Real}
    acc = XAccumulator(zero(T))
    @inbounds for i in eachindex(x)
        accumulate!(acc, x[i] * y[i])
    end
    return float(acc)
end

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

# manually optimized routine
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
