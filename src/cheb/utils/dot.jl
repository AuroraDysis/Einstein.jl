using Xsum

# TODO: benchmark this against the naive implementation
# The size of the vectors should be same
function dot_xsum(x::StridedVector{T}, y::StridedVector{T}) where {T<:Real}
    acc = XAccumulator{T}()
    @inbounds for i in eachindex(x)
        accumulate!(s, x[i] * y[i])
    end
    return T(acc)
end

function dot_kbn(A::StridedVector{T}, B::StridedVector{T}) where {T<:Real}
    @inbounds begin
        c = zero(T)
        s = A[1] * B[1] - c
        for i in 2:length(A)
            si = A[i] * B[i]
            t = s + si
            if abs(s) >= abs(si)
                c -= ((s - t) + si)
            else
                c -= ((si - t) + s)
            end
            s = t
        end
        return s - c
    end
end
