module PDESuiteTests

using Aqua
using JET
using PDESuite
using InlineTest

@testset "PDESuite.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(PDESuite)
    end

    @testset "Code linting (JET.jl)" begin
        JET.test_package(PDESuite; target_defined_modules=true)
    end
end

end

using ReTest
using PDESuite

retest(PDESuite, PDESuiteTests)
