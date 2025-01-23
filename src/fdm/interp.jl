# function fdm_interpwts(acc_order::Integer, dx::T) where {T<:AbstractFloat}
#     @argcheck acc_order > 0 "Accuracy order must be positive"
#     @argcheck dx > 0 "Spacing must be positive"

#     if acc_order == 1
#         return [0.0, 1.0]
#     end

#     num_pts = acc_order + 1
#     num_side = div(num_pts, 2)
#     local_grid = -num_side:num_side
#     wts = fdm_fornbergwts(0, 0, -n:0; hermite=false)

#     return wts
# end
