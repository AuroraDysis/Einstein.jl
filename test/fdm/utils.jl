@testitem "fdm_coeffnum" begin
    using GRSuite.FDMSuite, Test

    # Basic derivative orders with 2nd order accuracy
    @test fdm_coeffnum(1, 2) == 3
    @test fdm_coeffnum(2, 2) == 3
    @test fdm_coeffnum(3, 2) == 5

    # Higher derivatives with 4th order accuracy
    @test fdm_coeffnum(4, 4) == 7
    @test fdm_coeffnum(5, 4) == 9

    # Different accuracy orders for 1st derivative
    @test fdm_coeffnum(1, 2) == 3
    @test fdm_coeffnum(1, 4) == 5
    @test fdm_coeffnum(1, 6) == 7
    @test fdm_coeffnum(1, 8) == 9
end

@testitem "fdm_hermite_coeffnum" begin
    using GRSuite.FDMSuite, Test

    @test fdm_hermite_coeffnum(2, 4) == 3
    @test fdm_hermite_coeffnum(2, 8) == 5
end
