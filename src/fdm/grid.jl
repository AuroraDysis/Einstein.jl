function fdm_uniform_grid(lower_bound::TF, upper_bound::TF, dx::TF) where {TF<:AbstractFloat}
    @argcheck upper_bound > lower_bound "Invalid interval"
    @argcheck dx > 0 "Spacing must be positive"

    n = round(Int, (upper_bound - lower_bound) / dx) + 1
    precise_upper_bound = lower_bound + (n - 1) * dx

    # relative error check
    rtol = abs(precise_upper_bound - upper_bound) / abs(upper_bound)
    @argcheck rtol < typetol(TF) "Grid endpoint mismatch: |upper_bound - precise_upper_bound| = $(abs(upper_bound - precise_upper_bound)) exceeds tolerance ($(10 * eps(TF))). Consider adjusting dx to ensure upper_bound is reached precisely."

    return range(lower_bound, upper_bound, n)
end

export fdm_uniform_grid
