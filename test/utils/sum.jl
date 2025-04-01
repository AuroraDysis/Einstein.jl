@testitem "sum" begin
    function naive_sum(arr)
        sum = zero(eltype(arr))
        for i in eachindex(arr)
            sum += arr[i]
        end
        return sum
    end

    accurate_sum_functions = [sum_xsum, sum_kahan]

    for sum_func in accurate_sum_functions
        @testset "$sum_func" begin
            # 1. Summing many small numbers
            # Expected result: 1.0. Naive sum can lose precision as the sum grows.
            @testset "Summing many small numbers" begin
                n_small = 10_000_000
                small_val = 1.0e-7
                arr_small = fill(small_val, n_small)
                exact_small = Float64(n_small) * small_val

                @test sum_func(arr_small) ≈ exact_small
                @test naive_sum(arr_small) ≠ exact_small
            end

            # 2. Adding a small number to a large number repeatedly
            # The small number might be completely lost in naive summation.
            @testset "Adding a small number to a large number repeatedly" begin
                large_val = 1.0e18
                tiny_val = 1.0
                n_tiny = 1000
                arr_large_small = vcat([large_val], fill(tiny_val, n_tiny))
                # Exact sum is large_val + n_tiny * tiny_val
                exact_large_small = large_val + Float64(n_tiny) * tiny_val

                @test sum_func(arr_large_small) ≈ exact_large_small
                @test naive_sum(arr_large_small) ≠ exact_large_small
            end

            # 3. Ill-conditioned sum: Terms are large, s is small
            # Example: s = 1 + 1e16 - 1e16
            # Note: julia `sum` function is not able to handle this case well.
            @testset "Ill-conditioned sum" begin
                ill_val = 1.0e16
                arr_ill = [1.0, ill_val, -ill_val]
                # Exact sum is 1.0
                exact_ill = 1.0

                @test sum_func(arr_ill) ≈ exact_ill
                @test naive_sum(arr_ill) ≠ exact_ill
            end

            # 4. Normal case
            @testset "Normal case" begin
                arr_normal = [1.0, 2.0, 3.0]
                # Exact sum is 6.0
                exact_normal = 6.0

                @test sum_func(arr_normal) ≈ exact_normal
                @test naive_sum(arr_normal) ≈ exact_normal
            end
        end
    end
end
