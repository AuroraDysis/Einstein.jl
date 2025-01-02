# """
#     rectdiff2(m, n)

# Julia equivalent of the MATLAB function `rectdiff2(m, n)` for constructing
# a 1st-order rectangular differentiation matrix mapping from a 2nd-kind grid
# to a 1st-kind grid.

# Arguments:
# - `m` : Size of the output grid (number of rows).
# - `n` : Size of the input grid (number of columns).

# Returns:
# - `D` : The rectangular differentiation matrix (size `m x n`).
# """
# function rectdiff2(::Type{TR}, m::TI, n::TI) where {TR<:AbstractFloat,TI<:Integer}
#     nm1 = n - 1
#     cm1 = nm1 - m

#     t = cheb2_pts(n)    # Second-kind grid points
#     tau = cheb1_pts(m)  # First-kind grid points
#     T = cheb2_angles(n) # Second-kind angles
#     TAU = cheb1_angles(m) # First-kind angles

#     # Denominator using broadcasting
#     denom = [2 * sin((t + tau) / 2) * sin((tau - t) / 2) for tau in TAU, t in T]

#     # Numerator calculation
#     numer = [
#         (1 - tau * t) * (cos(cm1 * tau_angle) / sin(tau_angle)) for
#         (tau, tau_angle) in zip(tau, TAU), t in t
#     ]

#     sgn = (-1)^cm1

#     if cm1 == 0
#         D = numer ./ (denom .^ 2) / nm1
#     else
#         sin_part = sin.(cm1 .* TAU)
#         D = (sin_part .* ones(n)') ./ denom + numer ./ (denom .^ 2) / nm1
#         D .*= sgn
#     end

#     # Scale first and last columns
#     D[:, [1, n]] .*= 0.5

#     # if cm1 == 0
#     #     # Flipping trick
#     #     ii = rot180(tril(ones(Bool, m, n)))
#     #     rot_D = rot180(D)
#     #     D[ii] .= sgn .* rot_D[ii]
#     # end

#     # Sign adjustments
#     D[1:2:end, 1:2:end] .*= -1
#     D[2:2:end, 2:2:end] .*= -1

#     # Negative sum trick
#     # idx = [argmin(abs.(row)) for row in eachrow(denom)]
#     # for i in 1:m
#     #     D[i, idx[i]] = 0
#     #     D[i, idx[i]] = -sum(D[i, :])
#     # end

#     # if cm1 == 0
#     #     # Fix corner values
#     #     if m > 1 && nm1 != 0
#     #         corner_val = -0.25 / (nm1 * sin(π / (2 * m)) * sin(π / (4 * m))^2)
#     #         D[1, 1] = corner_val
#     #         D[end, end] = -corner_val

#     #         # Negative sum trick for corner entries
#     #         if n > 2
#     #             D[1, 2] = -sum(D[1, [1; 3:n]])
#     #             D[end, end - 1] = -D[1, 2]
#     #         end
#     #     end
#     # end

#     return D
# end

# export rectdiff2

# @testset "rectdiff2" begin
#     tol = 100 * eps(Float64)

#     D45 = rectdiff2(Float64, 4, 5)
#     D45_res = [
#         -5.5 6.82842712474619 -2.0 1.17157287525381 -0.5
#         -0.25 -1.35355339059327 2.0 -0.646446609406727 0.25
#         -0.25 0.646446609406727 -2.0 1.35355339059327 0.25
#         0.5 -1.17157287525381 2.0 -6.82842712474619 5.5
#     ]
#     println(D45 .- D45_res)
#     @test isapprox(D45, D45_res; atol=tol)

#     D56 = rectdiff2(Float64, 5, 6)
#     D56_res = [
#         -8.5 10.4721359549996 -2.89442719099992 1.52786404500042 -1.10557280900008 0.5
#         -0.865685424949239 -1.3023992240591 2.76251184000825 -0.931915350991663 0.603173584940985 -0.265685424949238
#         -0.1 0.305572809000084 -2.09442719099992 2.09442719099992 -0.305572809000084 0.1
#         0.265685424949238 -0.603173584940985 0.931915350991663 -2.76251184000825 1.3023992240591 0.865685424949239
#         -0.5 1.10557280900008 -1.52786404500042 2.89442719099992 -10.4721359549996 8.5
#     ]
#     println(D56 .- D56_res)
#     @test isapprox(D56, D56_res; atol=tol)
# end
