@testitem "dot" begin
    # Define a naive dot product function for comparison
    function naive_dot(a, b)
        # Ensure inputs are compatible
        if length(a) != length(b)
            throw(DimensionMismatch("vectors must have same length"))
        end
        # Use standard floating-point arithmetic
        s = zero(promote_type(eltype(a), eltype(b)))
        for i in eachindex(a)
            s += a[i] * b[i]
        end
        return s
    end

    # List the accurate dot product functions to test
    accurate_dot_functions = [dot_xsum, dot_kahan]

    for dot_func in accurate_dot_functions
        @testset "$dot_func" begin
            # 1. Dot product equivalent to summing many small numbers
            # We achieve this by dotting the array of small numbers with an array of ones.
            # Expected result: 1.0. Naive dot product (summation part) can lose precision.
            @testset "Dot product of many small numbers" begin
                n_small = 10_000_000
                small_val = 1.0e-7
                a_small = fill(small_val, n_small)
                b_small = ones(Float64, n_small) # Vector of ones
                exact_small = Float64(n_small) * small_val # Expected: 10_000_000 * 1.0e-7 = 1.0

                @test dot_func(a_small, b_small) ≈ exact_small
                @test naive_dot(a_small, b_small) != exact_small # Naive sum inside dot fails
            end

            # 2. Dot product equivalent to adding a small number to a large number repeatedly
            # Achieved by dotting an array [large, tiny, tiny,...] with an array of ones.
            # The small contributions might be lost in naive summation within the dot product.
            @testset "Dot product involving large and small numbers" begin
                large_val = 1.0e18
                tiny_val = 1.0
                n_tiny = 1000
                # Vector 'a' contains the mix of values
                a_large_small = vcat([large_val], fill(tiny_val, n_tiny))
                # Vector 'b' is ones to make dot(a,b) == sum(a)
                b_large_small = ones(Float64, length(a_large_small))
                # Exact sum is large_val + n_tiny * tiny_val
                exact_large_small = large_val + Float64(n_tiny) * tiny_val

                @test dot_func(a_large_small, b_large_small) ≈ exact_large_small
                # Naive sum inside dot is expected to fail here
                @test naive_dot(a_large_small, b_large_small) != exact_large_small
            end

            # 3. Ill-conditioned dot product: Sum involves large terms cancelling out
            # Example: dot([1.0, 1e16, -1e16], [1.0, 1.0, 1.0]) = 1*1 + 1e16*1 - 1e16*1 = 1.0
            # Naive summation within the dot product can fail.
            @testset "Ill-conditioned dot product" begin
                ill_val = 1.0e16
                # Vector 'a' has the ill-conditioned terms
                a_ill = [1.0, ill_val, -ill_val]
                # Vector 'b' is ones to make dot(a,b) == sum(a)
                b_ill = ones(Float64, length(a_ill))
                # Exact result is 1.0
                exact_ill = 1.0

                @test dot_func(a_ill, b_ill) ≈ exact_ill
                # Naive sum inside dot expected to fail (1.0 + 1e16 = 1e16, then 1e16 - 1e16 = 0.0)
                @test naive_dot(a_ill, b_ill) != exact_ill
            end

            # 4. Normal case dot product
            # Use vectors where the naive approach should also work.
            @testset "Normal case dot product" begin
                a_normal = [1.0, 2.0, 3.0]
                b_normal = [4.0, 5.0, 6.0]
                # Exact dot product = 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32.0
                exact_normal = 32.0

                @test dot_func(a_normal, b_normal) ≈ exact_normal
                @test naive_dot(a_normal, b_normal) ≈ exact_normal
            end
        end
    end
end