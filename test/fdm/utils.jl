@testitem "fdm_centralnum" begin
    using Einstein.FDMSuite, Test

    # Basic derivative orders with 2nd order accuracy
    @test fdm_centralnum(1, 2) == 3
    @test fdm_centralnum(2, 2) == 3
    @test fdm_centralnum(3, 2) == 5

    # Higher derivatives with 4th order accuracy
    @test fdm_centralnum(4, 4) == 7
    @test fdm_centralnum(5, 4) == 9

    # Different accuracy orders for 1st derivative
    @test fdm_centralnum(1, 2) == 3
    @test fdm_centralnum(1, 4) == 5
    @test fdm_centralnum(1, 6) == 7
    @test fdm_centralnum(1, 8) == 9
end

@testitem "fdm_hermitenum" begin
    using Einstein.FDMSuite, Test

    @test fdm_hermitenum(2, 4) == 3
    @test fdm_hermitenum(2, 8) == 5
end
