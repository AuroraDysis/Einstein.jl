module FiniteDifferenceSuite

using InlineTest

using ArgCheck: @argcheck

include("grid.jl")
include("weights.jl")
include("dissipation.jl")

end
