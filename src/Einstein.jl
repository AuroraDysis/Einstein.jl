module Einstein

using Reexport

include("utils/utils.jl")
@reexport using .Utils

include("cheb/cheb.jl")
@reexport using .ChebyshevSuite

include("fdm/fdm.jl")
@reexport using .FiniteDifferenceSuite

include("qnm/qnm.jl")
@reexport using .QNMSuite

using PrecompileTools

PrecompileTools.@compile_workload begin
    using .ChebyshevSuite

    cheb1_points(Float64, 5)
    cheb2_points(Float64, 5)

    cheb1_angles(Float64, 5)
    cheb2_angles(Float64, 5)

    # TODO: implement the rest of the precompiles
end

end
