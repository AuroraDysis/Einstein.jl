using TestItems

@testitem "GaussChebyshevGrid - angles" begin
    angles_0 = Float64[]
    angles_1 = [1.570796326794897]
    angles_2 = [2.356194490192345, 0.785398163397448]
    angles_5 = [2.827433388230814, 2.199114857512855, 1.570796326794897, 0.942477796076938, 0.314159265358979]

    @test GaussChebyshevGrid.angles(Float64, 0) ≈ angles_0
    @test GaussChebyshevGrid.angles(Float64, 1) ≈ angles_1
    @test GaussChebyshevGrid.angles(Float64, 2) ≈ angles_2
    @test GaussChebyshevGrid.angles(Float64, 5) ≈ angles_5
end
