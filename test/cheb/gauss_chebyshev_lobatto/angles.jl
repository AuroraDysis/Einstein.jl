using TestItems

@testitem "GaussChebyshevLobatto - angles" begin
    angles_0 = Float64[]
    angles_1 = [1.570796326794897]
    angles_2 = [3.141592653589793, 0.0]
    angles_5 = [
        3.141592653589793, 2.356194490192345, 1.570796326794897, 0.785398163397448, 0.0
    ]

    @test GaussChebyshevLobatto.angles(Float64, 0) ≈ angles_0
    @test GaussChebyshevLobatto.angles(Float64, 1) ≈ angles_1
    @test GaussChebyshevLobatto.angles(Float64, 2) ≈ angles_2
    @test GaussChebyshevLobatto.angles(Float64, 5) ≈ angles_5
end
