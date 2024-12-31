using GRSuite
using SafeTestsets
using Test
using Aqua
using JET

@testset "GRSuite.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(GRSuite)
    end

    @testset "Code linting (JET.jl)" begin
        JET.test_package(GRSuite; target_defined_modules=true)
    end

    @time @safetestset "Chebyshev pseudospectral method" include("cheb_test.jl")
    @time @safetestset "Finite difference method" include("fdm_test.jl")
end
