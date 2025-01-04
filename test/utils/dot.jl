@testitem "dot_xsum" begin
    using GRSuite.Utils, Test

    x = rand(10)
    y = rand(10)
    @test dot_xsum(x, y) ≈ sum(x .* y)

    x = [1.0, 2.0, 3.0]
    y = [4.0, 5.0, 6.0]
    @test dot_xsum(x, y) ≈ 32.0
end

@testitem "dot_kahan" begin
    using GRSuite.Utils, Test

    x = rand(10)
    y = rand(10)
    @test dot_kahan(x, y) ≈ sum(x .* y)

    x = [1.0, 2.0, 3.0]
    y = [4.0, 5.0, 6.0]
    @test dot_kahan(x, y) ≈ 32.0
end

@testitem "dot_kahan_neumaier" begin
    using GRSuite.Utils, Test

    x = rand(10)
    y = rand(10)
    @test dot_kahan_neumaier(x, y) ≈ sum(x .* y)

    x = [1.0, 2.0, 3.0]
    y = [4.0, 5.0, 6.0]
    @test dot_kahan_neumaier(x, y) ≈ 32.0
end
