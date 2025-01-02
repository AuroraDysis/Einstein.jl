# """
#     rectdiff1(m, n)

# Julia equivalent of the MATLAB function `rectdiff1(m, n)` for constructing
# a 1st-order rectangular differentiation matrix mapping from a 1st-kind grid.

# Arguments:
# - `m` : Size of the output grid (number of rows).
# - `n` : Size of the input grid (number of columns).

# Returns:
# - `D` : The rectangular differentiation matrix (size `m x n`).
# """
# function rectdiff1(::Type{TR}, m::TI, n::TI) where {TR<:AbstractFloat,TI<:Integer}
#     # mapping-from grid (angles):
#     T = cheb1_angles(n)'        # Row vector of length n
#     # difference between dimensions:
#     c = n - m
#     # mapping-to grid (angles):
#     TAU = cheb1_angles(m)       # Column vector of length m

#     # Denominator term
#     denom = @. 2 * sin((T + TAU) / 2) * sin((TAU - T) / 2)

#     # Sign matrix
#     sgn = ones(TR, m, n)
#     if isodd(c)
#         sgn[1:2:m, 2:2:n] .= -1
#         sgn[2:2:m, 1:2:n] .= -1
#     else
#         sgn[1:2:m, 1:2:n] .= -1
#         sgn[2:2:m, 2:2:n] .= -1
#     end

#     # Matrix construction (fixed broadcasting)
#     D = @. sgn * ((cos(c * TAU) / sin(TAU)) * sin(T) - sin(c * TAU) / n * sin(T) / denom) /
#         denom

#     # Negative-sum trick
#     # idx = [argmin(abs.(denom[i, :])) for i in 1:m]
#     # for i in 1:m
#     #     D[i, idx[i]] = 0
#     #     D[i, idx[i]] = -sum(@view D[i, :])
#     # end

#     # Flipping trick
#     # if c == 1
#     #     rot_D_180 = rot180(D)
#     #     mask = rot180(tril(ones(Bool, m, n)))
#     #     D[mask] .-= rot_D_180[mask]
#     # end

#     return D
# end

# export rectdiff1

# @testset "rectdiff1" begin
#     tol = 100 * eps(Float64)

#     D45 = rectdiff1(Float64, 4, 5)
#     D45_res = [
#         -4.57116486120251 6.35944256641146 -2.70279423558635 1.31914227171008 -0.404625741332680
#         -0.0484519638509037 -1.91971151095453 2.34412101225198 -0.504026090801922 0.128068553355378
#         -0.128068553355378 0.504026090801922 -2.34412101225198 1.91971151095453 0.0484519638509037
#         0.404625741332680 -1.31914227171008 2.70279423558635 -6.35944256641146 4.57116486120251
#     ]
#     @test isapprox(D45, D45_res; atol=tol)
#     println(D45 .- D45_res)

#     D56 = rectdiff1(Float64, 5, 6)
#     D56_res = [
#         -6.71910606480245 9.53284806144464 -4.39831637342411 2.49110926664733 -1.32569196231062 0.419157072445230
#         -0.253223738073512 -2.39105495582767 3.33681201184053 -1.01065674661971 0.453608214480045 -0.135484785799679
#         -0.0462335694140099 0.235702260395516 -2.40325617336917 2.40325617336917 -0.235702260395516 0.0462335694140099
#         0.135484785799679 -0.453608214480045 1.01065674661971 -3.33681201184053 2.39105495582767 0.253223738073512
#         -0.419157072445230 1.32569196231062 -2.49110926664733 4.39831637342411 -9.53284806144464 6.71910606480245
#     ]
#     @test isapprox(D56, D56_res; atol=tol)
#     println(D56 .- D56_res)
# end
