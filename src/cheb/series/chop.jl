"""
    ChebyshevSeriesChopContext{TF<:AbstractFloat, TI<:Integer}

A context object for the `cheb_series_chop!` function to reuse temporary storage
and minimize dynamic memory allocations.

Create instances using `cheb_series_chop_context`.
"""
struct ChebyshevSeriesChopContext{TF<:AbstractFloat,TI<:Integer}
    n_max::TI # Maximum coefficient vector length supported by this context
    abs_coeffs::Vector{TF}
    envelope::Vector{TF}
    log_env_vals::Vector{TF}
    objective_vals::Vector{TF}

    """
        ChebyshevSeriesChopContext{TF, TI}(n_max::TI) where {TF<:AbstractFloat, TI<:Integer}

    Internal constructor for the context. Use `cheb_series_chop_context`.
    """
    function ChebyshevSeriesChopContext{TF,TI}(
        n_max::TI
    ) where {TF<:AbstractFloat,TI<:Integer}
        # The algorithm itself requires n >= 17, but the context can be smaller
        # if only used for smaller inputs (though chopping won't occur then).
        # We enforce n_max >= 1 for valid vector indexing.
        @argcheck n_max >= 1 "Context maximum size n_max must be at least 1."
        abs_coeffs = Vector{TF}(undef, n_max)
        envelope = Vector{TF}(undef, n_max)
        log_env_vals = Vector{TF}(undef, n_max)
        objective_vals = Vector{TF}(undef, n_max)
        return new{TF,TI}(n_max, abs_coeffs, envelope, log_env_vals, objective_vals)
    end
end

"""
    cheb_series_chop_context(::Type{TF}, n_max::TI) where {TF<:AbstractFloat, TI<:Integer}

Create a context for `cheb_series_chop!` that can handle coefficient vectors
up to length `n_max`.

# Arguments
- `TF`: The floating-point type (e.g., `Float64`).
- `n_max`: The maximum length of the coefficient vector the context should support.

# Returns
- `ChebyshevSeriesChopContext`: The reusable context object.
"""
function cheb_series_chop_context(
    ::Type{TF}, n_max::TI
) where {TF<:AbstractFloat,TI<:Integer}
    return ChebyshevSeriesChopContext{TF,TI}(n_max)
end

function _cheb_series_chop_tol_impl!(
    coeffs::AbstractVector{TFC}, tol::TF
) where {TF<:AbstractFloat,TFC<:Union{TF,Complex{TF}}}
    @inbounds for k in length(coeffs):-1:1
        if abs(coeffs[k]) > tol
            return k
        end
    end
    return 0
end

"""
    _cheb_series_chop_standard_impl!(ctx::ChebyshevSeriesChopContext{TF}, coeffs::AbstractVector{TF}, tol::TF) where {TF<:AbstractFloat}

Internal core logic for `cheb_series_chop!`, operating using the pre-allocated `ctx`.
Returns the cutoff index. This function performs the actual computation.
"""
function _cheb_series_chop_standard_impl!(
    ctx::ChebyshevSeriesChopContext{TF,TI}, coeffs::AbstractVector{TF}, tol::TF
) where {TF<:AbstractFloat,TI<:Integer}
    n = length(coeffs)
    (; abs_coeffs, envelope, log_env_vals, objective_vals) = ctx

    # Other checks (coeffs empty, tol range) are done in the calling functions.

    # --- Use views into pre-allocated vectors for the current size n ---
    # These views allow operating on the relevant portion of the context's arrays
    # without allocating new memory for each call.
    b_view = view(abs_coeffs, 1:n)         # For abs(coeffs)
    envelope_view = view(envelope, 1:n)    # For envelope calculation and storage
    log_env_view = view(log_env_vals, 1:n) # For log10(envelope) storage (up to search_end)
    obj_vals_view = view(objective_vals, 1:n) # For objective function storage (up to search_end)

    # Initialize cutoff to n (default "not happy" state, meaning no chop point found)
    cutoff = n

    # --- Step 1: Compute envelope (using views and in-place operations) ---
    # Calculate absolute values of coefficients and store in b_view
    @. b_view = abs(coeffs)

    # Compute reverse cumulative maximum into envelope_view.
    # envelope_view[j] will store max(abs(coeffs[k])) for k >= j.
    if n > 0 # Ensure there's at least one element
        envelope_view[n] = b_view[n] # Last element is its own max
        # Iterate backwards to compute cumulative maximum efficiently
        for j in (n - 1):-1:1
            envelope_view[j] = max(b_view[j], envelope_view[j + 1])
        end
    end

    # Normalize the envelope by its first element (the maximum absolute coefficient).
    # Store the normalized envelope back into envelope_view.
    m1 = envelope_view[1]
    if m1 <= 0.0 # Use <= to handle potential -0.0
        # If the largest coefficient is zero (or negative zero), all are zero.
        return 1 # Chop immediately after the first coefficient.
    end
    # Perform in-place normalization.
    @. envelope_view = envelope_view / m1

    # --- Step 2: Scan for plateau in the normalized envelope ---
    # A plateau indicates where the coefficients stop decaying significantly.
    plateau_found::Bool = false
    plateau_point::Int = 0 # Index where the plateau starts
    plateau_j2::Int = 0    # Index used for the plateau check calculation
    log_tol::TF = log(tol) # Precompute log(tol) for efficiency

    # Iterate through potential start points of the plateau
    for j in 2:n
        # Calculate the lookahead index based on the heuristic from the paper/chebfun
        current_j2 = round(Int, 1.25 * j + 5)
        # Stop if the lookahead index goes beyond the array bounds
        if current_j2 > n
            break
        end

        # Get envelope values at the current point and the lookahead point
        e1 = envelope_view[j]
        e2 = envelope_view[current_j2]

        # Check the plateau condition
        local plateau::Bool # Ensure 'plateau' has local scope within the loop iteration
        if e1 == 0.0
            # If envelope value is zero, we've definitely hit a plateau (or the end of non-zeros)
            plateau = true
        else
            # Heuristic condition from the reference algorithm
            r::TF = 3.0 * (1.0 - log(e1) / log_tol)
            plateau = (e2 / e1 > r) # Compare ratio to the threshold 'r'
        end

        # If a plateau is detected, record the indices and stop scanning
        if plateau
            plateau_point = j - 1 # Plateau starts just before the current index 'j'
            plateau_j2 = current_j2
            plateau_found = true
            break # Exit the loop once the first plateau is found
        end
    end

    # If no plateau was found after scanning, we can't chop.
    if !plateau_found
        return n # Return original length, indicating no chop point found
    end

    # --- Step 3: Refine cutoff point based on the plateau ---
    # Ensure the plateau point is valid (should always be positive if found)
    if plateau_point <= 0
        # This case should ideally not happen if plateau_found is true.
        # Add a warning and return a safe default.
        @warn "Internal logic error: Plateau point is non-positive ($plateau_point) despite plateau being found. Returning cutoff=1."
        return 1
    end

    # If the envelope value at the plateau point is exactly zero,
    # it means all subsequent coefficients are also zero (due to cumulative max).
    if envelope_view[plateau_point] == 0.0
        cutoff = plateau_point # Chop right at the start of the zero region
    else
        # If the envelope is non-zero, we need to find the minimum of an objective function
        # to balance accuracy and series length.

        # Define a target value based on the tolerance
        target_val::TF = tol^(7.0 / 6.0) # Heuristic target value

        # Find the last index j3 where the envelope value is >= target_val.
        # This helps limit the search range for the minimum.
        j3 = findlast(x -> x >= target_val, view(envelope_view, 1:plateau_point)) # Search only up to plateau_point
        if isnothing(j3)
            # If no value meets the target, search starts from the beginning.
            # However, the logic below uses plateau_j2 as the primary upper bound.
            # We refine the search end based on j3 relative to plateau_j2.
            j3 = 0 # Set to 0 if not found
        end

        # Determine the end index for the minimization search range.
        # Start with the lookahead index from the plateau scan.
        search_end::Int = plateau_j2
        adjust_objective_end::Bool = false
        # If j3 is found and is earlier than plateau_j2, we can potentially shorten the search range.
        if j3 > 0 && j3 < plateau_j2
            # Search up to j3 + 1 to include the transition point.
            search_end = j3 + 1
            # Flag that the last value of the objective function might need adjustment
            # because we truncated the search based on target_val.
            adjust_objective_end = true
        end
        # Ensure search_end is within valid bounds [1, n]
        search_end = max(1, search_end) # Must be at least 1
        search_end = min(search_end, n)   # Cannot exceed total length

        # If the search range is empty (e.g., if plateau_j2 was < 1 somehow), chop at 1.
        if search_end <= 0 # Should technically be search_end < 1, but <= 0 is safer
            cutoff = 1
        else
            # --- Minimize the objective function using pre-allocated views ---
            # Create views for the active search range [1:search_end]
            current_env_active_view = view(envelope_view, 1:search_end)
            log_env_active_view = view(log_env_view, 1:search_end) # View for log10(env)
            obj_vals_active_view = view(obj_vals_view, 1:search_end) # View for objective func

            # Calculate log10 of the envelope in-place into log_env_active_view.
            # Handle cases where envelope value is 0.
            @. log_env_active_view = ifelse(
                current_env_active_view <= 0.0, -TF(Inf), log10(current_env_active_view)
            )

            # Calculate the bias term for the objective function.
            # LinRange allocation is relatively minor compared to the main loop, kept for clarity.
            log10_tol::TF = log10(tol)
            # Bias decreases linearly from 0 to -log10(tol)/3 over the search range.
            bias = LinRange(TF(0.0), (-TF(1.0) / TF(3.0)) * log10_tol, search_end)

            # Calculate the objective function: log10(envelope) + bias. Store in-place.
            @. obj_vals_active_view = log_env_active_view + bias

            # If we truncated the search range based on j3, adjust the last objective value.
            # This prevents artificially choosing the end point just because we stopped searching early.
            if adjust_objective_end && search_end > 0
                # Calculate the log10 of the target value used for truncation.
                log_adjusted_val::TF = (target_val <= 0.0) ? -TF(Inf) : log10(target_val)
                # Replace the last objective value with the value corresponding to target_val.
                # Use the precomputed bias value at the end of the (truncated) range.
                obj_vals_active_view[end] = log_adjusted_val + bias[end]
            end

            # Find the index of the minimum value in the objective function view.
            # Manual loop is often slightly more performant than findmin for simple cases
            # and avoids allocating a tuple for (value, index).
            min_val::TF = obj_vals_active_view[1]
            min_idx::Int = 1
            # Use @inbounds because search_end is validated against array bounds.
            @inbounds for i in 2:search_end
                current_val = obj_vals_active_view[i]
                if current_val < min_val
                    min_val = current_val
                    min_idx = i
                end
            end

            # The optimal cutoff index is one less than the index of the minimum objective value.
            # Ensure the cutoff is at least 1.
            cutoff = max(min_idx - 1, 1)
        end
    end

    # Return the calculated cutoff index.
    return cutoff
end

"""
    cheb_series_chop!(ctx::ChebyshevSeriesChopContext{TF}, coeffs::AbstractVector{TF}, tol::TF=eps(TF)) where {TF<:AbstractFloat}

Determine the chopping index using a pre-allocated context.

# Arguments
- `ctx`: The `ChebyshevSeriesChopContext` object.
- `coeffs`: The coefficient vector (length must be `<= ctx.n_max`).
- `tol`: The relative tolerance (must be `0 < tol < 1`). Defaults to `eps(TF)`.

# Returns
- `cutoff::Int`: The calculated cutoff index. If `cutoff == length(coeffs)`, no suitable chopping point was found.
"""
function cheb_series_chop!(
    ctx::ChebyshevSeriesChopContext{TF,TI}, coeffs::AbstractVector{TF}, tol::TF=eps(TF)
) where {TF<:AbstractFloat,TI<:Integer}
    @boundscheck begin
        @argcheck !isempty(coeffs) "coeffs must not be empty"
        @argcheck 0 < tol < 1 "tol must be between 0 and 1 (exclusive)"
        @argcheck length(coeffs) <= ctx.n_max "coeffs must be <= ctx.n_max"
    end

    n = length(coeffs)

    # Algorithm requires a minimum length to work reliably.
    # If the input is too short, we cannot reliably find a chop point.
    if n < 17
        return _cheb_series_chop_tol_impl!(coeffs, tol)
    end

    return _cheb_series_chop_standard_impl!(ctx, coeffs, tol)
end

"""
    cheb_series_chop(coeffs::AbstractVector{TF}, tol::TF = eps(TF)) where {TF<:AbstractFloat}

Determine a suitable cutoff index for a coefficient vector using the "standard" chopping rule [aurentz2015choppingchebyshevseries](@cite).

This version creates a temporary context internally and is suitable for single use or when performance is not critical.
For repeated use with vectors of similar maximum size, create and reuse a
`ChebyshevSeriesChopContext` via `cheb_series_chop_context` and call `cheb_series_chop!` for better performance (reduced allocations).

# Arguments
- `coeffs::AbstractVector{TF}`: A non-empty vector of coefficients.
- `tol::TF`: A target relative accuracy, a number in `(0, 1)`. Defaults to machine epsilon (`eps(TF)`).

# Returns
- `cutoff::Int`: A positive integer representing the chopping point.
    - If `cutoff == length(coeffs)`, a satisfactory chopping point was not found.
    - If `cutoff < length(coeffs)`, this is the recommended last index to keep (i.e., `coeffs[1:cutoff]`).

# Examples
```julia
using ArgCheck # Required for the function if not loaded globally

coeffs = 10.0 .^ -(0:49) # Example decaying coefficients
random_noise = cos.((1.0:50.0).^2) # Example noise

# --- Using the convenience function (allocates context internally) ---
println("--- Convenience Function Examples ---")
cutoff1 = cheb_series_chop(coeffs) # Default tolerance
println("Example 1 (coeffs only): ", cutoff1) # Expected: ~18

cutoff2 = cheb_series_chop(coeffs .+ 1e-16 .* random_noise)
println("Example 2 (coeffs + 1e-16 noise): ", cutoff2) # Expected: ~15

cutoff3 = cheb_series_chop(coeffs .+ 1e-13 .* random_noise)
println("Example 3 (coeffs + 1e-13 noise): ", cutoff3) # Expected: ~13

cutoff4 = cheb_series_chop(coeffs .+ 1e-10 .* random_noise)
println("Example 4 (coeffs + 1e-10 noise): ", cutoff4) # Expected: 50 (noise dominates)

cutoff5 = cheb_series_chop(coeffs .+ 1e-10 .* random_noise, 1e-10) # Specify tolerance
println("Example 5 (coeffs + 1e-10 noise, tol=1e-10): ", cutoff5) # Expected: ~10

# --- Using a pre-allocated context (more performant for repeated calls) ---
println("\\n--- Pre-allocated Context Examples ---")
max_len = 60 # Max size the context should handle
ctx = cheb_series_chop_context(Float64, max_len)

coeffs_a = 10.0 .^ -(0:49)
coeffs_b = 10.0 .^ -(0:39) .+ 1e-14 .* cos.((1.0:40.0).^2)

# Use the context via the cheb_series_chop! function
cutoff_ctx1 = cheb_series_chop!(ctx, coeffs_a) # Default tolerance
println("Context Example 1 (coeffs_a): ", cutoff_ctx1) # Expected: ~18

cutoff_ctx2 = cheb_series_chop!(ctx, coeffs_b, 1e-12) # Specify tolerance
println("Context Example 2 (coeffs_b, tol=1e-12): ", cutoff_ctx2) # Expected: ~14
```

# Performance Guide
For best performance, especially in loops or repeated calls:
```julia
# 1. Create context once
ctx = cheb_series_chop_context(Float64, n_max) # Choose appropriate TF and n_max

# 2. Reuse context in loops
for coeffs_vec in list_of_coefficient_vectors
    cutoff = cheb_series_chop!(ctx, coeffs_vec)
    # Use the cutoff...
end
```

# References
- [aurentz2015choppingchebyshevseries](@citet*) Aurentz, J. L., & Trefethen, L. N. (2015). Chopping a Chebyshev series. ACM Transactions on Mathematical Software (TOMS), 41(4), 1-18.
- [chebfun/standardChop.m](https://github.com/chebfun/chebfun/blob/master/%40chebtech/standardChop.m) (Implementation in Chebfun, MATLAB)
"""
function cheb_series_chop(
    coeffs::AbstractVector{TF}, tol::TF=eps(TF)
) where {TF<:AbstractFloat}
    @boundscheck begin
        @argcheck !isempty(coeffs) "coeffs must not be empty"
        @argcheck 0 < tol < 1 "tol must be between 0 and 1 (exclusive)"
    end

    n = length(coeffs)

    # Algorithm requires a minimum length to work reliably.
    # If the input is too short, we cannot reliably find a chop point.
    if n < 17
        return _cheb_series_chop_tol_impl!(coeffs, tol)
    end

    ctx = cheb_series_chop_context(TF, n)

    return _cheb_series_chop_standard_impl!(ctx, coeffs, tol)
end

export cheb_series_chop, cheb_series_chop!, cheb_series_chop_context

#= Example Usage (Illustrative - run in a Julia environment)
using ArgCheck # Make sure ArgCheck is loaded

# Define coefficients and noise as in the docstring example
coeffs = 10.0 .^ -(0:49)
random_noise = cos.((1.0:50.0).^2)

# --- Convenience Function ---
println("--- Convenience Function Examples ---")
println("Example 1: ", cheb_series_chop(coeffs))
println("Example 2: ", cheb_series_chop(coeffs .+ 1e-16 .* random_noise))
println("Example 3: ", cheb_series_chop(coeffs .+ 1e-13 .* random_noise))
println("Example 4: ", cheb_series_chop(coeffs .+ 1e-10 .* random_noise))
println("Example 5: ", cheb_series_chop(coeffs .+ 1e-10 .* random_noise, 1e-10))

# --- Pre-allocated Context ---
println("\n--- Pre-allocated Context Examples ---")
max_len = 60
ctx = cheb_series_chop_context(Float64, max_len)

coeffs_a = 10.0 .^ -(0:49)
coeffs_b = 10.0 .^ -(0:39) .+ 1e-14 .* cos.((1.0:40.0).^2)

cutoff_ctx1 = cheb_series_chop!(ctx, coeffs_a)
println("Context Example 1 (coeffs_a): ", cutoff_ctx1)

cutoff_ctx2 = cheb_series_chop!(ctx, coeffs_b, 1e-12)
println("Context Example 2 (coeffs_b, tol=1e-12): ", cutoff_ctx2)

# --- Benchmarking (Optional - requires BenchmarkTools) ---
# using BenchmarkTools
# println("\n--- Allocation Benchmarks ---")
# print("Convenience function: ")
# @btime cheb_series_chop($coeffs)
# print("Context function:     ")
# @btime cheb_series_chop!($ctx, $coeffs)

=#