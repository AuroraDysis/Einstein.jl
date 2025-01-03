using TestItems
using TestItemRunner

@testitem "Code quality (Aqua.jl)" begin
    using Aqua, GRSuite

    Aqua.test_all(GRSuite)
end

@testitem "Code linting (JET.jl)" begin
    using JET, GRSuite

    JET.test_package(GRSuite; target_defined_modules=true)
end

@run_package_tests
