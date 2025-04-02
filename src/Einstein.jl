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
    for grid in (ChebyshevGaussGrid, ChebyshevLobattoGrid)
        points = grid.points(Float64, 5)
        angles = grid.angles(Float64, 5)
    end

    # TODO: implement the rest of the precompiles
end

end
