module Einstein

using Reexport

include("utils/utils.jl")
@reexport using .Utils

include("cheb/cheb.jl")
@reexport using .ChebSuite

include("fdm/fdm.jl")
@reexport using .FDMSuite

include("qnm/qnm.jl")
@reexport using .QNMSuite

using PrecompileTools

PrecompileTools.@compile_workload begin
    using .ChebSuite

    cheb1_pts(Float64, 5)
    cheb2_pts(Float64, 5)

    cheb1_angles(Float64, 5)
    cheb2_angles(Float64, 5)

    # TODO: implement the rest of the precompiles
end

end
