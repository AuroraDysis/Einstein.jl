using GRSuite
using Test
using Aqua
using JET

@testset "GRSuite.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(GRSuite)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(GRSuite; target_defined_modules = true)
    end
    # Write your tests here.
end
