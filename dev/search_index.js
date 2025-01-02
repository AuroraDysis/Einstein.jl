var documenterSearchIndex = {"docs":
[{"location":"cheb/#Chebyshev-Suite","page":"Chebyshev Suite","title":"Chebyshev Suite","text":"","category":"section"},{"location":"cheb/","page":"Chebyshev Suite","title":"Chebyshev Suite","text":"Modules = [PDESuite.ChebyshevSuite]","category":"page"},{"location":"cheb/","page":"Chebyshev Suite","title":"Chebyshev Suite","text":"Modules = [PDESuite.ChebyshevSuite]","category":"page"},{"location":"cheb/#PDESuite.ChebyshevSuite.Cheb1Coeffs2ValsOp","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.Cheb1Coeffs2ValsOp","text":"cheb1_coeffs2vals(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}\nop::Cheb1Coeffs2ValsOp([TR=Float64], n::TI)(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR},TI<:Integer}\n\nConvert Chebyshev coefficients to values at Chebyshev points of the 1st kind.\n\nPerformance Guide\n\nFor best performance, especially in loops or repeated calls:\n\nop = Cheb1Coeffs2ValsOp(Float64, n)\nvalues = op(coeffs)\n\nReferences\n\nchebfun/@chebtech1/coeffs2vals.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"type"},{"location":"cheb/#PDESuite.ChebyshevSuite.Cheb1Vals2CoeffsOp","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.Cheb1Vals2CoeffsOp","text":"cheb1_vals2coeffs(vals::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}\nop::Cheb1Vals2CoeffsOp([TR=Float64], n::TI)(vals::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR},TI<:Integer}\n\nConvert values at Chebyshev points of the 1st kind into Chebyshev coefficients.\n\nPerformance Guide\n\nFor best performance, especially in loops or repeated calls:\n\nop = Cheb1Vals2CoeffsOp(Float64, n)\nvalues = op(coeffs)\n\nReferences\n\nchebfun/@chebtech1/vals2coeffs.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"type"},{"location":"cheb/#PDESuite.ChebyshevSuite.Cheb2Vals2CoeffsOp","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.Cheb2Vals2CoeffsOp","text":"cheb2_vals2coeffs(vals::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}\nop::Cheb2Vals2CoeffsOp([TR=Float64], n::TI)(vals::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}\n\nConvert values at Chebyshev points of the 2nd kind into Chebyshev coefficients.\n\nPerformance Guide\n\nFor best performance, especially in loops or repeated calls:\n\nop = Cheb2Vals2CoeffsOp(Float64, n)\nvalues = op(coeffs)\n\nReferences\n\nchebfun/@chebtech2/vals2coeffs.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"type"},{"location":"cheb/#PDESuite.ChebyshevSuite.ChebCumsumOp","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.ChebCumsumOp","text":"cheb_cumsum(f::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}\nChebCumsumOp([TR=Float64], n::TI)(f::VT) where {TR<:AbstractFloat,TI<:Integer,VT<:AbstractVector{TR}}\n\nCompute the indefinite integral of a function given its Chebyshev coefficients.\n\nArguments\n\nf: Vector of Chebyshev coefficients of the function to be integrated\n\nReferences\n\nchebfun/@chebtech/cumsum.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"type"},{"location":"cheb/#PDESuite.ChebyshevSuite.bary-Union{Tuple{VT3}, Tuple{VT2}, Tuple{VT1}, Tuple{TR}, Tuple{VT1, VT2, VT3, TR}} where {TR<:AbstractFloat, VT1<:AbstractVector{TR}, VT2<:AbstractVector{TR}, VT3<:AbstractVector{TR}}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.bary","text":"bary(w::VT1, x::VT2, f::VT3, x0::TR) where {\n    TR<:AbstractFloat,\n    VT1<:AbstractVector{TR},\n    VT2<:AbstractVector{TR},\n    VT3<:AbstractVector{TR},\n}\n\nEvaluate a polynomial interpolant using the barycentric interpolation formula.\n\nArguments\n\nw: Vector of barycentric weights\nx: Vector of interpolation points (typically Chebyshev points)\nf: Vector of function values at interpolation points\nx0: Point at which to evaluate the interpolant\n\nReference\n\nchebfun/@chebtech2/bary.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.bary_diffmat-Union{Tuple{TI}, Tuple{VT3}, Tuple{VT2}, Tuple{VT1}, Tuple{TR}, Tuple{VT1, VT2, TI, VT3}} where {TR<:AbstractFloat, VT1<:AbstractVector{TR}, VT2<:AbstractVector{TR}, VT3<:AbstractVector{TR}, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.bary_diffmat","text":"bary_diffmat(x; w=nothing, k=1, t=nothing)\n\nCompute the barycentric differentiation matrix.\n\nReferences:\n\nchebfun/@chebcolloc/baryDiffMat.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb1_amat-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb1_amat","text":"cheb1_amat([TR=Float64], n::Integer)\n\nConstruct the analysis matrix A that transforms function values at Chebyshev points of the 1st kind to Chebyshev coefficients.\n\nArguments\n\nTR: Element type (defaults to Float64)\nn: Number of points/coefficients\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb1_angles-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb1_angles","text":"cheb1_angles([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}\n\nCompute angles for Chebyshev points of the 1st kind: theta_k = frac(2k + 1)pi2n quad k = n-1ldots0\n\nArguments\n\nTR: Type parameter for the angles (e.g., Float64)\nn: Number of points\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb1_barywts-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb1_barywts","text":"cheb1_barywts([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}\n\nCompute the barycentric weights for Chebyshev points of the 1st kind.\n\nArguments\n\nTR: Type parameter for the weights (e.g., Float64)\nn: Number of points\n\nReferences\n\nBerrut and Trefethen [BT04]\nchebfun/@chebtech1/barywts.m at master · chebfun/chebfun\n\nSee also: bary, cheb1_pts\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb1_cumsummat-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb1_cumsummat","text":"cheb1_cumsummat([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}\ncheb1_cumsummat([TR=Float64], n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}\n\nCompute Chebyshev integration matrix that maps function values at n Chebyshev points of the 1st kind to values of the integral of the interpolating polynomial at those points, with the convention that the first value is zero.\n\nReferences\n\nchebfun/@chebcolloc1/chebcolloc1.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb1_diffmat-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}, Tuple{Type{TR}, TI, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb1_diffmat","text":"cheb1_diffmat([TR=Float64], n::TI, k::TI=1) where {TR<:AbstractFloat,TI<:Integer}\n\nConstruct a Chebyshev differentiation that maps function values at n Chebyshev points of the 1st kind  to values of the k-th derivative of the interpolating polynomial at those points.\n\nArguments\n\nTR: Element type (defaults to Float64)\nn::Integer: Number of Chebyshev points\nk::Integer=1: Order of the derivative (default: 1)\n\nReferences\n\nchebfun/@chebcolloc1/chebcolloc1.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb1_pts-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb1_pts","text":"cheb1_pts([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}\ncheb1_pts([TR=Float64], n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}\n\nGenerate Chebyshev points of the 2nd kind.\n\nFor the standard interval [-1,1]: x_k = -cosleft(frac(2k + 1)pi2nright) quad k = 01ldotsn-1\n\nFor mapped interval [xmin,xmax]: x_mapped = fracx_max + x_min2 + fracx_min - x_max2x_k\n\nArguments\n\nTR: Type parameter for the grid points (e.g., Float64)\nn: Number of points\nx_min: (Optional) Lower bound of the mapped interval\nx_max: (Optional) Upper bound of the mapped interval\n\nReferences\n\nchebfun/@chebtech1/chebpts.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb1_quadwts-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb1_quadwts","text":"cheb1_quadwts([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}\n\nCompute quadrature weights for Chebyshev points of the 1st kind.\n\nArguments\n\nTR: Type parameter for the weights (e.g., Float64)\nn: Number of points\n\nReferences\n\nchebfun/@chebtech1/quadwts.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb1_smat-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb1_smat","text":"cheb1_smat([TR=Float64], n::Integer)\n\nConstruct the synthesis matrix S that transforms Chebyshev coefficients to function values at Chebyshev points of the 1st kind.\n\nArguments\n\nTR: Element type (defaults to Float64)\nn: Number of points/coefficients\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb2_amat-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb2_amat","text":"cheb2_amat([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}\n\nConstruct the analysis matrix A that transforms function values at Chebyshev points of the 2nd kind to Chebyshev coefficients.\n\nArguments\n\nTR: Element type (defaults to Float64)\nn: Number of points/coefficients\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb2_angles-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb2_angles","text":"cheb2_angles([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}\n\nCompute angles for Chebyshev points of the 2nd kind: theta_k = frackpin-1 quad k = n-1ldots0\n\nArguments\n\nTR: Type parameter for the angles (e.g., Float64)\nn: Number of points\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb2_barywts-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb2_barywts","text":"cheb2_barywts([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}\n\nCompute the barycentric weights for Chebyshev points of the 2nd kind.\n\nArguments\n\nTR: Type parameter for the weights (e.g., Float64)\nn: Number of points\n\nReferences\n\nchebfun/@chebtech2/barywts.m at master · chebfun/chebfun\n\nSee also: bary, cheb2_pts\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb2_coeffs2vals-Union{Tuple{VT}, Tuple{TR}} where {TR<:AbstractFloat, VT<:AbstractVector{TR}}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb2_coeffs2vals","text":"cheb2_coeffs2vals(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}\nop::Cheb2Coeffs2ValsOp([TR=Float64], n::TI)(coeffs::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR},TI<:Integer}\n\nConvert Chebyshev coefficients to values at Chebyshev points of the 2nd kind.\n\nPerformance Guide\n\nFor best performance, especially in loops or repeated calls:\n\nop = Cheb2Coeffs2ValsOp(Float64, n)\nvalues = op(coeffs)\n\nReferences\n\nchebfun/@chebtech2/coeffs2vals.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb2_cumsummat-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb2_cumsummat","text":"cheb2_cumsummat([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}\ncheb2_cumsummat([TR=Float64], n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}\n\nCompute Chebyshev integration matrix that maps function values at n Chebyshev points of the 2st kind to values of the integral of the interpolating polynomial at those points, with the convention that the first value is zero.\n\nReferences\n\nchebfun/@chebcolloc2/chebcolloc2.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb2_diffmat-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}, Tuple{Type{TR}, TI, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb2_diffmat","text":"cheb2_diffmat([TR=Float64], n::TI, k::TI=1) where {TR<:AbstractFloat,TI<:Integer}\n\nConstruct a Chebyshev differentiation that maps function values at n Chebyshev points of the 2nd kind  to values of the k-th derivative of the interpolating polynomial at those points.\n\nArguments\n\nTR: Element type (defaults to Float64)\nn::Integer: Number of Chebyshev points\nk::Integer=1: Order of the derivative (default: 1)\n\nReferences\n\nchebfun/@chebcolloc2/chebcolloc2.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb2_pts-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb2_pts","text":"cheb2_pts([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}\ncheb2_pts([TR=Float64], n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}\n\nGenerate Chebyshev points of the 1st kind.\n\nFor the standard interval [-1,1]: x_k = -cosleft(frackpin-1right) quad k = 01ldotsn-1\n\nFor mapped interval [xmin,xmax]: x_mapped = fracx_max + x_min2 + fracx_min - x_max2x_k\n\nArguments\n\nTR: Type parameter for the grid points (e.g., Float64)\nn: Number of points\nx_min: (Optional) Lower bound of the mapped interval\nx_max: (Optional) Upper bound of the mapped interval\n\nReferences\n\nchebfun/@chebtech2/chebpts.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb2_quadwts-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb2_quadwts","text":"cheb2_quadwts([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}\n\nCompute quadrature weights for Chebyshev points of the 2nd kind.\n\nArguments\n\nTR: Type parameter for the weights (e.g., Float64)\nn: Number of points\n\nReferences\n\nchebfun/@chebtech2/quadwts.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb2_smat-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb2_smat","text":"cheb1_smat([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}\n\nConstruct the synthesis matrix S that transforms Chebyshev coefficients to function values at Chebyshev points of the 2nd kind.\n\nArguments\n\nTR: Element type (defaults to Float64)\nn: Number of points/coefficients\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb_clenshaw-Union{Tuple{VT}, Tuple{T}, Tuple{VT, T}} where {T<:AbstractFloat, VT<:AbstractVector{T}}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb_clenshaw","text":"cheb_clenshaw(c::VT, x::T) where {T<:AbstractFloat,VT<:AbstractVector{T}}\n\nEvaluate Chebyshev coefficients at a point using Clenshaw's algorithm.\n\nArguments\n\nc: Vector of Chebyshev coefficients c_0 c_1 ldots c_n\nx: Evaluation point in [-1,1]\n\nReferences\n\nchebfun/@chebtech/clenshaw.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb_coeffs_cumsummat-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb_coeffs_cumsummat","text":"cheb_coeffs_cumsummat([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}\ncheb_coeffs_cumsummat([TR=Float64], n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}\n\nGenerate the Chebyshev coefficient integration matrix that maps Chebyshev coefficients to the coefficients of the integral of the interpolating polynomial.\n\nArguments\n\nTR: Type parameter for the matrix elements (e.g., Float64)\nn: Size of the matrix (n×n)\nx_min: (Optional) Lower bound of the integration interval\nx_max: (Optional) Upper bound of the integration interval\n\nReferences\n\nchebfun/@chebcolloc1/chebcolloc1.m at master · chebfun/chebfun\nchebfun/@chebcolloc2/chebcolloc2.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb_feval-Union{Tuple{VT}, Tuple{TR}, Tuple{VT, TR}} where {TR<:AbstractFloat, VT<:AbstractVector{TR}}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb_feval","text":"cheb_feval(f::VT, x::TR) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}\n\nEvaluate Chebyshev coefficients at a point.\n\nPerformance Notes\n\nClenshaw's algorithm: O(n) operations per point\n(TODO) NDCT: O(n log n) operations for many points simultaneously\n\nReferences\n\nchebfun/@chebtech/feval.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb_rectdiff1-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb_rectdiff1","text":"cheb_rectdiff1([TR=Float64], m::TI, n::TI) where {TR<:AbstractFloat,TI<:Integer}\ncheb_rectdiff1([TR=Float64], m::TI, n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}\n\nConstructing a 1st-order rectangular differentiation matrix mapping from a 1st-kind grid\n\nArguments:\n\nm : Size of the output grid (number of rows).\nn : Size of the input grid (number of columns).\n\nReferences\n\nchebfun/diffmat.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb_rectdiff2-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb_rectdiff2","text":"cheb_rectdiff2([TR=Float64], m::TI, n::TI) where {TR<:AbstractFloat,TI<:Integer}\ncheb_rectdiff2([TR=Float64], m::TI, n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}\n\nConstruct a 1st-order rectangular differentiation matrix mapping from a 2nd-kind grid.\n\nArguments:\n\nm : Size of the output grid (number of rows)\nn : Size of the input grid (number of columns)\n\nReferences\n\nchebfun/diffmat.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.cheb_rectint-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.cheb_rectint","text":"cheb_rectint([TR=Float64], n::TI) where {TR<:AbstractFloat,TI<:Integer}\ncheb_rectint([TR=Float64], n::TI, x_min::TR, x_max::TR) where {TR<:AbstractFloat,TI<:Integer}\n\nGenerate the Chebyshev integration matrix that operates directly on function values.\n\nArguments\n\nTR: Type parameter for the matrix elements (e.g., Float64)\nn: Size of the matrix (n×n)\nx_min: (Optional) Lower bound of the integration interval\nx_max: (Optional) Upper bound of the integration interval\n\nReferences\n\nchebfun/intmat.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.runtests-Tuple","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.runtests","text":"PDESuite.ChebyshevSuite.runtests(pattern...; kwargs...)\n\nEquivalent to ReTest.retest(PDESuite.ChebyshevSuite, pattern...; kwargs...). This function is defined automatically in any module containing a @testset, possibly nested within submodules.\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.ultra_convertmat-Union{Tuple{TI}, Tuple{TR}, Tuple{Type{TR}, TI, TI, TI}} where {TR<:AbstractFloat, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.ultra_convertmat","text":"ultra_convertmat([TR=Float64], n::TI, K1::TI, K2::TI) where {TR<:AbstractFloat,TI<:Integer}\n\nConversion matrix used in the ultraspherical spectral method. Returns N-by-N matrix realization of conversion operator between ultraspherical polynomial bases. Maps N coefficients from C^{(K1)} basis to C^{(K2)} basis.\n\nReferences\n\nchebfun/@ultraS/convertmat.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.ultra_diffmat-Union{Tuple{TI}, Tuple{TI, TI}} where TI<:Integer","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.ultra_diffmat","text":"ultra_diffmat(n::TI, m::TI) where {TI<:Integer}\n\nDifferentiation matrices for ultraspherical spectral method that takes n Chebyshev coefficients and returns n C^(m) coefficients that represent the derivative of the Chebyshev series. Here, C^(k) is the ultraspherical polynomial basis with parameter k.\n\nArguments\n\nn::TI: Number of points\nm::TI: Order of differentiation\n\nReferences\n\nchebfun/@ultraS/diffmat.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.ultra_multmat-Union{Tuple{TI}, Tuple{VT}, Tuple{TR}, Tuple{VT, TI}} where {TR<:AbstractFloat, VT<:AbstractVector{TR}, TI<:Integer}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.ultra_multmat","text":"ultra_multmat(a::VT, λ::TI) where {TR<:AbstractFloat,VT<:AbstractVector{TR},TI<:Integer}\n\nConstruct nxn multiplication matrix representing the multiplication of F in the C^(lambda) basis.\n\nReferences\n\nchebfun/@ultraS/multmat.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.ultra_spconvert-Union{Tuple{TR}, Tuple{TI}, Tuple{Type{TR}, TI, TI}} where {TI<:Integer, TR<:AbstractFloat}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.ultra_spconvert","text":"ultra_spconvert([TR=Float64], n::TI, λ::TR) where {TR<:AbstractFloat,TI<:Integer,TR<:Real}\n\nCompute sparse representation for conversion operators. Returns the truncation of the operator that transforms C^lambda (Ultraspherical polynomials) to C^lambda + 1. The truncation gives back a matrix of size n x n.\n\nReferences\n\nchebfun/@ultraS/spconvert.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/#PDESuite.ChebyshevSuite.ultra_sphankel-Union{Tuple{VT}, Tuple{TR}} where {TR<:AbstractFloat, VT<:AbstractVector{TR}}","page":"Chebyshev Suite","title":"PDESuite.ChebyshevSuite.ultra_sphankel","text":"sphankel(r::VT) where {TR<:AbstractFloat,VT<:AbstractVector{TR}}\n\nConstruct a sparse Hankel matrix by forming it as an upside-down Toeplitz matrix. This is required by the ultraspherical multiplication operator.\n\nReferences\n\nchebfun/@ultraS/sphankel.m at master · chebfun/chebfun\n\n\n\n\n\n","category":"method"},{"location":"cheb/","page":"Chebyshev Suite","title":"Chebyshev Suite","text":"M. C. Babiuc and others. Implementation of standard testbeds for numerical relativity. Class. Quant. Grav. 25, 125012 (2008), arXiv:0709.3559 [gr-qc].\n\n\n\nJ.-P. Berrut and L. N. Trefethen. Barycentric lagrange interpolation. SIAM review 46, 501–517 (2004).\n\n\n\nB. Fornberg. Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of computation 51, 699–706 (1988).\n\n\n\nB. Fornberg. Classroom Note:Calculation of Weights in Finite Difference Formulas. SIAM Review 40, 685–691 (1998), arXiv:https://doi.org/10.1137/S0036144596322507.\n\n\n\nB. Fornberg. An algorithm for calculating Hermite-based finite difference weights. IMA Journal of Numerical Analysis 41, 801–813 (2021).\n\n\n\n","category":"page"},{"location":"fdm/#Finite-Difference-Suite","page":"Finite Difference Suite","title":"Finite Difference Suite","text":"","category":"section"},{"location":"fdm/","page":"Finite Difference Suite","title":"Finite Difference Suite","text":"Modules = [PDESuite.FiniteDifferenceSuite]","category":"page"},{"location":"fdm/","page":"Finite Difference Suite","title":"Finite Difference Suite","text":"Modules = [PDESuite.FiniteDifferenceSuite]","category":"page"},{"location":"fdm/#PDESuite.FiniteDifferenceSuite.dissipation_order-Tuple{TI} where TI<:Integer","page":"Finite Difference Suite","title":"PDESuite.FiniteDifferenceSuite.dissipation_order","text":"dissipation_order(acc_order::Integer)\n\nCalculate the order of dissipation needed for a given finite difference accuracy order [B+08]. For a scheme of accuracy order 2r-2, returns dissipation order 2r.\n\n\n\n\n\n","category":"method"},{"location":"fdm/#PDESuite.FiniteDifferenceSuite.dissipation_wts-Tuple{TI} where TI<:Integer","page":"Finite Difference Suite","title":"PDESuite.FiniteDifferenceSuite.dissipation_wts","text":"calculate_dissipation_weights(diss_order::Integer)\n\nCalculate the weights for Kreiss-Oliger dissipation of given order [B+08].\n\n\n\n\n\n","category":"method"},{"location":"fdm/#PDESuite.FiniteDifferenceSuite.fdm_grid-Union{Tuple{TR}, Tuple{Type{TR}, TR, TR, TR}} where TR<:AbstractFloat","page":"Finite Difference Suite","title":"PDESuite.FiniteDifferenceSuite.fdm_grid","text":"fdm_grid([TR=Float64], x_min::TR, x_max::TR, dx::TR) where {TR<:AbstractFloat}\n\nGenerate a uniform grid for finite difference methods.\n\nArguments\n\nTR: Type parameter for the grid points (e.g., Float64)\nx_min: Lower bound of the interval\nx_max: Upper bound of the interval\ndx: Grid spacing\n\nReturns\n\nVector of n uniformly spaced points, where n = round((xmax - xmin)/dx) + 1\n\n\n\n\n\n","category":"method"},{"location":"fdm/#PDESuite.FiniteDifferenceSuite.fornberg_calculate_wts-Union{Tuple{TI}, Tuple{VT}, Tuple{T2}, Tuple{T}, Tuple{TI, T, VT}} where {T<:Real, T2<:Real, VT<:AbstractVector{T2}, TI<:Integer}","page":"Finite Difference Suite","title":"PDESuite.FiniteDifferenceSuite.fornberg_calculate_wts","text":"fornberg_calculate_wts([T=Float64], order::Integer, x0::Real, x::AbstractVector; \n                         dfdx::Bool=false)\n\nCalculate finite difference weights for arbitrary-order derivatives using the Fornberg algorithm.\n\nArguments\n\nT: Type parameter for the weights (defaults to type of x0)\norder: Order of the derivative to approximate\nx0: Point at which to approximate the derivative\nx: Grid points to use in the approximation\ndfdx: Whether to include first derivative values (Hermite finite differences)\n\nReturns\n\nIf dfdx == false:\n\nVector{T}: Weights for function values\n\nIf dfdx == true:\n\nTuple{Vector{T}, Vector{T}}: Weights for (function values, derivative values)\n\nMathematical Background\n\nFor a function f(x), the derivative approximation takes the form:\n\nIf dfdx == false (standard finite differences):\n\nf^(n)(x_0) approx sum_j=1^N c_j f(x_j)\n\nIf dfdx == true (Hermite finite differences):\n\nf^(n)(x_0) approx sum_j=1^N d_j f(x_j) + e_j f(x_j)\n\nRequirements\n\nFor standard finite differences: N > order\nFor Hermite finite differences: N > order/2 + 1\n\nwhere N is the length of x\n\nExamples\n\n# Standard central difference for first derivative\nx = [-1.0, 0.0, 1.0]\nw = fornberg_calculate_wts(1, 0.0, x)\n# Returns approximately [-0.5, 0.0, 0.5]\n\n# Forward difference for second derivative\nx = [0.0, 1.0, 2.0, 3.0]\nw = fornberg_calculate_wts(2, 0.0, x)\n\n# Hermite finite difference for third derivative\nx = [-1.0, 0.0, 1.0]\nw_f, w_d = fornberg_calculate_wts(3, 0.0, x, dfdx=true)\n\nReferences\n\nFornberg [For88]\nFornberg [For21]\nFornberg [For98]\nMethodOfLines.jl/src/discretization/schemes/fornbergcalculatewts.jl at master · SciML/MethodOfLines.jl\nprecision - Numerical derivative and finite difference coefficients: any update of the Fornberg method? - Computational Science Stack Exchange\n\nNotes\n\nThe implementation includes a stability correction for higher-order derivatives\nFor first derivatives (order=1), the weights sum to zero\nThe algorithm handles both uniform and non-uniform grids\nWhen using Hermite finite differences, fewer points are needed but derivatives must be available\n\nSee also: fdm_grid\n\n\n\n\n\n","category":"method"},{"location":"fdm/#PDESuite.FiniteDifferenceSuite.runtests-Tuple","page":"Finite Difference Suite","title":"PDESuite.FiniteDifferenceSuite.runtests","text":"PDESuite.FiniteDifferenceSuite.runtests(pattern...; kwargs...)\n\nEquivalent to ReTest.retest(PDESuite.FiniteDifferenceSuite, pattern...; kwargs...). This function is defined automatically in any module containing a @testset, possibly nested within submodules.\n\n\n\n\n\n","category":"method"},{"location":"fdm/","page":"Finite Difference Suite","title":"Finite Difference Suite","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = PDESuite","category":"page"},{"location":"#PDESuite","page":"Home","title":"PDESuite","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for PDESuite.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [PDESuite]","category":"page"},{"location":"#PDESuite.runtests-Tuple","page":"Home","title":"PDESuite.runtests","text":"PDESuite.runtests(pattern...; kwargs...)\n\nEquivalent to ReTest.retest(PDESuite, pattern...; kwargs...). This function is defined automatically in any module containing a @testset, possibly nested within submodules.\n\n\n\n\n\n","category":"method"}]
}
