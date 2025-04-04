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
    n = 3
    for grid in (ChebyshevGaussGrid, ChebyshevLobattoGrid)
        points = grid.points(n)
        angles = grid.angles(n)
        coeffs = grid.coeffs2vals(points)
        vals = grid.vals2coeffs(coeffs)
        plan_coeffs2vals = grid.coeffs2vals_plan(n)
        plan_vals2coeffs = grid.vals2coeffs_plan(n)
        coeffs_from_plan = plan_coeffs2vals * coeffs
        vals_from_plan = plan_vals2coeffs * vals
        coeffs_matrix = grid.coeffs2vals_matrix(n)
        vals_matrix = grid.vals2coeffs_matrix(n)
        weights = grid.barycentric_weights(n)
        quadrature_weights = grid.quadrature_weights(n)
        differentiation_matrix = grid.differentiation_matrix(n)
        integration_matrix = grid.integration_matrix(n)
    end

    # TODO: implement the rest of the precompiles
end

end
