using ReTest, PDESuite

include("PDESuiteTests.jl")

retest(PDESuite, PDESuiteTests)
