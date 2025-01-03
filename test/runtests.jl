using Test
using Aqua
using JET
using PDESuite
using SafeTestsets

@testset "PDESuite.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(PDESuite)
    end

    @testset "Code linting (JET.jl)" begin
        JET.test_package(PDESuite; target_defined_modules=true)
    end

    @safetestset "ChebSuite" begin
        include("cheb/cheb.jl")
    end
    @safetestset "FDMSuite" begin
        include("fdm/fdm.jl")
    end
end
