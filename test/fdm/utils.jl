@testitem "fdm_central_width" begin
    # Basic derivative orders with 2nd order accuracy
    @test fdm_central_width(1, 2) == 3
    @test fdm_central_width(2, 2) == 3
    @test fdm_central_width(3, 2) == 5

    # Higher derivatives with 4th order accuracy
    @test fdm_central_width(4, 4) == 7
    @test fdm_central_width(5, 4) == 9

    # Different accuracy orders for 1st derivative
    @test fdm_central_width(1, 2) == 3
    @test fdm_central_width(1, 4) == 5
    @test fdm_central_width(1, 6) == 7
    @test fdm_central_width(1, 8) == 9
end

@testitem "fdm_hermite_width" begin
    @test fdm_hermite_width(2, 4) == 3
    @test fdm_hermite_width(2, 8) == 5
end
