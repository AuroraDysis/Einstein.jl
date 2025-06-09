using TestItems
using TestItemRunner

@testitem "Code quality (Aqua.jl)" begin
    using Aqua

    Aqua.test_all(Einstein)
end

@testitem "Code linting (JET.jl)" begin
    using JET

    JET.test_package(Einstein; target_defined_modules=true)
end

@run_package_tests
