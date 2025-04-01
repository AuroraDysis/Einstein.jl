using TestItems

@testitem "GaussChebyshevLobatto - barycentric_weights" begin
    weights_0 = Float64[]
    weights_1 = [1.0]
    weights_2 = [-0.500000000000000, 0.500000000000000]
    weights_5 = [0.500000000000000, -1.0, 1.0, -1.0, 0.500000000000000]
    weights_6 = [-0.500000000000000, 1.0, -1.0, 1.0, -1.0, 0.500000000000000]

    @test GaussChebyshevLobatto.barycentric_weights(0) ≈ weights_0
    @test GaussChebyshevLobatto.barycentric_weights(1) ≈ weights_1
    @test GaussChebyshevLobatto.barycentric_weights(2) ≈ weights_2
    @test GaussChebyshevLobatto.barycentric_weights(5) ≈ weights_5
    @test GaussChebyshevLobatto.barycentric_weights(6) ≈ weights_6
end
