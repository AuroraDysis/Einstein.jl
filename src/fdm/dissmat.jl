"""
    fdm_dissmat(::Type{TR}, diss_order::Integer, n::Integer; transpose::Bool=false)

Create a finite difference dissipation matrix with specified dissipation order.

# Arguments
- `TR`: The element type of the matrix
- `diss_order::Integer`: The order of the dissipation operator
- `n::Integer`: The size of the matrix
- `transpose::Bool=false`: Whether to return the transpose of the matrix
"""
function fdm_dissmat(
    ::Type{TR}, diss_order::Integer, n::Integer; transpose::Bool=false
) where {TR<:Real}
    op = fdm_dissop(diss_order, one(TR), one(TR))
    num_side = op.num_side
    wts = op.wts

    dissmat = BandedMatrix(Zeros{TR}(n, n), (num_side, num_side))

    @inbounds for i in (num_side + 1):(n - num_side)
        dissmat[i, (i - num_side):(i + num_side)] .= wts
    end

    if transpose
        dissmat_T = similar(dissmat)
        transpose!(dissmat_T, dissmat)
        return dissmat_T
    else
        return dissmat
    end
end

export fdm_dissmat
