using TestItems
using TestItemRunner

@testitem "Code quality (Aqua.jl)" begin
    using Aqua, PDESuite

    Aqua.test_all(PDESuite)
end

@testitem "Code linting (JET.jl)" begin
    using JET, PDESuite

    JET.test_package(PDESuite; target_defined_modules=true)
end

@run_package_tests
