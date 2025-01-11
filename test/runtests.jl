using TestItems
using TestItemRunner

@testitem "Code quality (Aqua.jl)" begin
    using Aqua, Einstein

    Aqua.test_all(Einstein)
end

@run_package_tests
