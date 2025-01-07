using TestItems
using TestItemRunner

@testitem "Code quality (Aqua.jl)" begin
    using Aqua, GRSuite

    Aqua.test_all(GRSuite)
end

@run_package_tests
