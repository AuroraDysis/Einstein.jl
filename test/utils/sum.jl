@testitem "sum_xsum" begin
    using GRSuite.Utils, Test

    v = rand(10)
    @test sum_xsum(v) ≈ sum(v)

    v = [1.0, 2.0, 3.0]
    @test sum_xsum(v) ≈ 6.0
end

@testitem "sum_kahan" begin
    using GRSuite.Utils, Test

    v = rand(10)
    @test sum_kahan(v) ≈ sum(v)

    v = [1.0, 2.0, 3.0]
    @test sum_kahan(v) ≈ 6.0
end

@testitem "sum_kahan_opt" begin
    using GRSuite.Utils, Test

    v = rand(10)
    @test sum_kahan_opt(v) ≈ sum(v)

    v = [1.0, 2.0, 3.0]
    @test sum_kahan_opt(v) ≈ 6.0
end