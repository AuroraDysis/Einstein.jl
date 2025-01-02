function cheb2_diffmat(::Type{TR}, n::TI, k::TI=1) where {TR<:AbstractFloat,TI<:Integer}
    x = cheb2_pts(TR, n)               # First kind points.
    w = cheb2_barywts(TR, n)           # Barycentric weights.
    t = cheb2_angles(TR, n)            # acos(x).
    D = bary_diffmat(x, w, k, t)       # Construct matrix.
    return D
end
