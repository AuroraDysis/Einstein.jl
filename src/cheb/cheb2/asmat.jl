"""
    cheb2_amat([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}

Construct the analysis matrix A that transforms function values at Chebyshev points of the 2nd kind to Chebyshev coefficients.

# Arguments
- `TR`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb2_amat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    A = Array{Float64,2}(undef, n, n)
    op = Cheb2Vals2CoeffsOp(TR, n)
    @inbounds for i in 1:n
        A[:, i] = op(OneElement(one(TR), i, n))
    end
    return A
end

function cheb2_amat(n::TI) where {TI<:Integer}
    return cheb2_amat(Float64, n)
end

"""
    cheb1_smat([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}

Construct the synthesis matrix S that transforms Chebyshev coefficients to function values at Chebyshev points of the 2nd kind.

# Arguments
- `TR`: Element type (defaults to Float64)
- `n`: Number of points/coefficients
"""
function cheb2_smat(::Type{TR}, n::TI) where {TR<:AbstractFloat,TI<:Integer}
    S = Array{Float64,2}(undef, n, n)
    op = Cheb2Coeffs2ValsOp(TR, n)
    @inbounds for i in 1:n
        S[:, i] = op(OneElement(one(TR), i, n))
    end
    return S
end

function cheb2_smat(n::TI) where {TI<:Integer}
    return cheb2_smat(Float64, n)
end

export cheb2_amat, cheb2_smat
