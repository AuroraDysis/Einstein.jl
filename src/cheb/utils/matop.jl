struct MatrixOp{TR<:AbstractFloat}
    m::Integer
    n::Integer
    rows::Vector{Vector{TR}}

    function MatrixOp(matrix::AbstractMatrix{TR}) where {TR<:AbstractFloat}
        m = size(matrix, 1)
        n = size(matrix, 2)
        rows = Vector{Vector{TR}}(undef, m)
        @inbounds for i in 1:m
            rows[i] = matrix[i, :]
        end
        return new{TR}(m, n, rows)
    end
end

function (op::MatrixOp)(row::Integer, vec::StridedVector{TR}) where {TR<:IEE}
    return dot_xsum(op.rows[row], vec)
end
