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

@testitem "sum_kahan_neumaier" begin
    using GRSuite.Utils, Test

    v = rand(10)
    @test sum_kahan_neumaier(v) ≈ sum(v)

    v = [1.0, 2.0, 3.0]
    @test sum_kahan_neumaier(v) ≈ 6.0
end