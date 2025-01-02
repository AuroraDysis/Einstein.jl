using ToeplitzMatrices

"""
    sphankel(r::VT) where {VT<:AbstractVector{TR},TR<:AbstractFloat}

this forms a sparse Hankel matrix by forming it as an upside-
down Toeplitz matrix. This is required by the ultraspherical multiplication
operator.
"""
function ultra_sphankel(r::VT) where {VT<:AbstractVector{TR},TR<:AbstractFloat}
    return Hankel(r, OneElement(r[end], 1, length(r)))
end

export ultra_sphankel
