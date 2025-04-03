"""
    cheb_series_chop(coeffs::AbstractVector{TF}, tol::TF = eps(TF)) where {TF<:AbstractFloat}

Determine a suitable cutoff index for a coefficient vector using the "standard" chopping rule [aurentz2015choppingchebyshevseries](@cite).

# Arguments
- `coeffs::AbstractVector{TF}`: A non-empty vector of coefficients (real or complex).
- `tol::TF`: A target relative accuracy, a number in `(0, 1)`. Defaults to machine epsilon (`eps(Float64)`).

# Returns
- `cutoff::Int`: A positive integer representing the chopping point.
    - If `cutoff == length(coeffs)`, a satisfactory chopping point was not found (the algorithm is "not happy").
    - If `cutoff < length(coeffs)`, this is the recommended last index to keep (the algorithm is "happy").

# Examples
```julia
coeffs = 10.0 .^ -(1:50)
random = cos.((1.0:50.0).^2)
println("Example 1: ", cheb_series_chop(coeffs)) # Should be 18
println("Example 2: ", cheb_series_chop(coeffs .+ 1e-16 .* random)) # Should be 15
println("Example 3: ", cheb_series_chop(coeffs .+ 1e-13 .* random)) # Should be 13
println("Example 4: ", cheb_series_chop(coeffs .+ 1e-10 .* random)) # Should be 50
println("Example 5: ", cheb_series_chop(coeffs .+ 1e-10 .* random, 1e-10)) # Should be 10
```

# References
- [aurentz2015choppingchebyshevseries](@citet*)
- [chebfun/standardChop.m at master · chebfun/chebfun](https://github.com/chebfun/chebfun/blob/master/standardChop.m)
"""
function cheb_series_chop(
    coeffs::AbstractVector{TF}, tol::TF=eps(TF)
) where {TF<:AbstractFloat}
    @boundscheck begin
        @argcheck !isempty(coeffs) "coeffs must not be empty"
        @argcheck 0 < tol < 1 "tol must be between 0 and 1"
    end

    n = length(coeffs)
    # Initialize cutoff to n (default "not happy" state)
    cutoff = n

    # Algorithm requires a minimum length to work reliably
    if n < 17
        return cutoff # Return n, indicating no chopping point found
    end

    # --- Step 1: Compute the monotonically nonincreasing envelope ---
    # Calculate absolute values of coefficients
    b = abs.(coeffs)

    # Compute the reverse cumulative maximum envelope
    # Initialize envelope vector with the appropriate float type
    m = Vector{TF}(undef, n)
    if n > 0 # Ensure coeffs is not empty (though n < 17 check handles this)
        m[n] = TF(b[n]) # Start from the last element
        # Iterate backwards to compute cumulative maximum
        for j in (n - 1):-1:1
            m[j] = max(TF(b[j]), m[j + 1])
        end
    end

    # Normalize the envelope by its first element
    m1 = m[1]
    if m1 == 0.0
        # If the maximum absolute coefficient is zero, all are zero. Chop at 1.
        return 1
    end
    envelope::Vector{TF} = m ./ m1 # Type assertion for clarity

    # --- Step 2: Scan the envelope for a 'plateau' ---
    # A plateau indicates a region where coefficients stop decaying significantly.
    plateau_found::Bool = false
    plateau_point::Int = 0      # Index j-1 where plateau starts
    plateau_j2::Int = 0         # Index round(1.25*j + 5) used in plateau check
    log_tol::TF = log(TF(tol))  # Precompute log(tol) for efficiency

    # Iterate through potential plateau starting points
    for j in 2:n
        # Calculate the end index for the potential plateau check
        current_j2 = round(Int, 1.25 * j + 5)

        # If the end index is beyond the vector length, no further plateaus possible
        if current_j2 > n
            break
        end

        e1 = envelope[j]    # Envelope value at the start of the potential plateau
        e2 = envelope[current_j2] # Envelope value at the end

        # Determine if a plateau exists based on relative decay
        local plateau::Bool # Scope variable locally to the loop iteration
        if e1 == 0.0
            # If envelope hits zero, consider it a plateau
            plateau = true
        else
            # Calculate the flatness tolerance `r`
            # r depends on how close e1 is to the tolerance `tol`
            # r = 3 * (1 - log(e1)/log(tol)) ranges from 0 (e1=tol) to 3 (e1=1)
            r::TF = 3.0 * (1.0 - log(e1) / log_tol)
            # Check if the decay is flatter than required by r
            plateau = (e2 / e1 > r)
        end

        # If a plateau is detected, record its position and exit the scan
        if plateau
            plateau_point = j - 1      # Plateau starts at index j-1
            plateau_j2 = current_j2  # Store the j2 corresponding to this plateau
            plateau_found = true
            break
        end
    end

    # If no plateau was found after scanning the entire envelope, return n
    if !plateau_found
        return n
    end

    # --- Step 3: Refine the cutoff point ---
    # A plateau was found starting at `plateau_point`. Now refine the exact `cutoff`.

    # Safety check (should be >= 1 given loop starts at j=2)
    if plateau_point <= 0
        @warn "Internal error: Plateau point is non-positive ($plateau_point). Returning 1."
        return 1
    end

    # If the envelope at the plateau start is already zero, chop there
    if envelope[plateau_point] == 0.0
        cutoff = plateau_point
    else
        # Define a threshold based on tolerance for refining the search range
        target_val::TF = TF(tol)^(7.0 / 6.0)

        # Find j3: the index of the last coefficient envelope value >= target_val
        j3 = findlast(x -> x >= target_val, envelope)
        # If no element meets the criteria, findlast returns nothing
        if isnothing(j3)
            j3 = 0
        end

        # Determine the end index (`search_end`) for the minimization search range [1...search_end]
        search_end::Int = plateau_j2 # Start with the j2 from Step 2
        adjust_objective_end::Bool = false # Flag if we need to adjust the objective function later
        if j3 < plateau_j2
            # If j3 is smaller, limit the search range based on j3
            search_end = j3 + 1
            adjust_objective_end = true
        end
        # Ensure search_end is within valid bounds [1, n]
        search_end = max(1, search_end)
        search_end = min(search_end, n)

        if search_end == 0 # Should not happen if n>=17
            cutoff = 1 # Fallback
        else
            # Minimize: log10(envelope) + linear bias term over the range [1...search_end]
            indices = 1:search_end
            # Use a view to avoid allocating a new vector for the slice
            env_vals_view = view(envelope, indices)

            # Calculate log10 of envelope values, handling log10(0) = -Inf
            log_env_vals = map(x -> x == 0.0 ? -Inf : log10(x), env_vals_view)

            # Calculate the linear bias term to penalize longer vectors
            log10_tol::TF = log10(TF(tol))
            # LinRange is inclusive of start and stop
            bias = LinRange(0.0, (-1.0 / 3.0) * log10_tol, search_end)

            # Combine log envelope and bias to get the objective function values
            objective_vals = log_env_vals .+ bias

            # Replicate MATLAB's adjustment: if search range was limited (j3 < plateau_j2),
            # it effectively sets envelope[search_end] = target_val before the log10/bias calculation.
            # We apply this adjustment directly to the last element of the objective function.
            if adjust_objective_end && search_end > 0
                # Value corresponding to envelope[search_end] = target_val
                log_adjusted_val::TF = log10(target_val)
                # Get the bias value for the last element
                # bias_last = bias[end] # LinRange ensures this is correct
                # Recalculate the last objective value
                # Handle log10(0) case for target_val if tol is extremely small
                objective_vals[end] =
                    (log_adjusted_val == -Inf ? -Inf : log_adjusted_val) + bias[end]
            end

            # Find the index of the minimum value in the objective function
            min_val, min_idx = findmin(objective_vals)

            # The cutoff is one less than the index of the minimum, but at least 1
            cutoff = max(min_idx - 1, 1)
        end
    end

    # Return the calculated cutoff index
    return cutoff
end

# # To make the function available if this code is part of a module:
# # export cheb_series_chop
# ```julia
# # Example Usage (as provided in the docstring)
# coeffs = 10.0 .^ -(1:50)
# random = cos.((1.0:50.0).^2)
# println("Example 1: ", cheb_series_chop(coeffs)) # Should be 18
# println("Example 2: ", cheb_series_chop(coeffs .+ 1e-16 .* random)) # Should be 15
# println("Example 3: ", cheb_series_chop(coeffs .+ 1e-13 .* random)) # Should be 13
# println("Example 4: ", cheb_series_chop(coeffs .+ 1e-10 .* random)) # Should be 50
# println("Example 5: ", cheb_series_chop(coeffs .+ 1e-10 .* random, 1e-10)) # Should be 10
# ```
