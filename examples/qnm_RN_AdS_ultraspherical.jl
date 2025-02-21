### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ fdf0edd0-efbc-11ef-1dea-e32569d65705
using ApproxFun, LinearAlgebra, PlutoUI, LaTeXStrings

# ╔═╡ 64dcba7b-daea-488c-bc89-e729cc7537e7
wrong(text) = Markdown.MD(Markdown.Admonition("danger", "Input is out of range!", [text]));

# ╔═╡ e57eb498-14bc-4606-89d0-3ddbaa20c5b5
md"""
We use `Float64` by default to obtain high accuracy results. For even higher precision, `Double64` (requires the `DoubleFloats.jl` package), `BigFloat`, or other types can be used.
"""

# ╔═╡ 6327b564-6066-45b4-b188-095e3761b96d
TF = Float64

# ╔═╡ b1db4296-da97-4abe-9303-b77878f212ac
md"""
# Calculate scalar QNMs of RN AdS black hole using the ultraspherical spectral method

We are interested in scalar perturbations of a Reissner-Nordström-AdS (RNAdS) black hole metric:

$$ds^2 = f(r)dt^2 - \frac{dr^2}{f(r)} - r^2(d\theta^2 + \sin^2\theta d\phi^2),$$ 

where $f(r) = \Delta/r^2$ and $\Delta = r^2 - 2Mr + Q^2 + r^4/R^2$. Here, $M$ denotes the black hole mass, $Q$ the charge, and $R$ the AdS radius.

## References

- [[gr-qc/0301052] Quasinormal modes of Reissner-Nordström-anti-de Sitter black holes: scalar, electromagnetic and gravitational perturbations](https://arxiv.org/abs/gr-qc/0301052)
- [[2209.09324] Calculating quasinormal modes of Schwarzschild anti-de Sitter black holes using the continued fraction method](https://arxiv.org/abs/2209.09324)
- [[hep-th/0003295] Quasinormal modes of Reissner-Nordström Anti-de Sitter Black Holes](https://arxiv.org/abs/hep-th/0003295)
- [[1202.1347] A fast and well-conditioned spectral method](https://arxiv.org/abs/1202.1347)
"""

# ╔═╡ 8d2ee7b2-e5e9-402a-95b0-b5c53db84c67
md"""
We choose units such that $R = 1$.
"""

# ╔═╡ 4d0cd86e-8cb4-4c40-a202-6513b6257892
R = one(TF)

# ╔═╡ 24060d90-60f9-4dd6-8206-a28514fc2ab1
md"""
**Input horizon radius $r_+$:**
"""

# ╔═╡ 63dc53d0-0c1c-432a-9a0e-07502660eb6d
@bind r₊_str TextField(default="1.0")

# ╔═╡ a5869853-874c-4779-866f-ae5e7af251a8
# parse input as a numeric value
r₊ = parse(TF, r₊_str)

# ╔═╡ 053e3de3-8ea4-409c-aa00-5d2dc07d8f14
md"""
The extremal value of the black hole charge, $Q_{\text{ext}}$, is given by the following function of the hole's horizon radius:

$$Q_{\text{ext}}^2 = r_{+}^2 \left(1 + \frac{3 r_{+}^2}{R^2}\right) \,.$$
"""

# ╔═╡ 07200c40-621c-46f0-ae95-58a3b4755648
begin
	Qext = r₊ * sqrt(1 + 3 * r₊^2 / R^2)
	L"""Q_\text{ext} = %$(Qext)"""
end

# ╔═╡ 33b34c18-883c-46e2-bc84-80f5e411d382
md"""
**Input black hole charge $Q$:**
"""

# ╔═╡ e687ba54-efd4-4d0d-967c-eb20160a8c51
@bind Q_str TextField(default="0.0")

# ╔═╡ b766f294-5c18-40ce-8f92-696c1dcaaacb
# parse input as a numeric value
begin
	Q = parse(TF, Q_str)

	if Q < Qext
		L"Q / Q_\text{ext} = %$(Q / Qext)"
	else
		wrong(md"Out of range! $Q$ should be greater than zero and smaller than $Q_\text{ext}$.")
	end
end

# ╔═╡ 7f6a4b13-c252-47e9-b01c-5db1a7710985
md"""
The black hole mass relates to its charge $Q$ and horizon radius $r_+$ via:  

$$M = \frac{1}{2}\left(r_+ + \frac{r_+^3}{R^2} + \frac{Q^2}{r_+}\right).$$
"""

# ╔═╡ a5dc2c36-0c26-4cfe-a03a-137aabb47467
M = (r₊ + r₊^3 / R^2 + Q^2 / r₊) / 2

# ╔═╡ 9a09cfdf-98ae-4c73-943a-91472c0798e8
md"""
The radial part of scalar perturbations satisfies the wave equation:  

$$\frac{d^2 \psi}{dr_*^2} + \left(\omega^2 - V\right)\psi = 0,$$  

with potential:

$$V = f(r)\left[\frac{l(l+1)}{r^2} + \frac{f'(r)}{r}\right].$$  

Boundary conditions enforce ingoing waves at the horizon and no outgoing waves at infinity:

$$\psi(r) \approx \begin{cases}  
e^{-i\omega r_*} & \text{as } r \to r_+, \\  
0 & \text{as } r \to \infty.  
\end{cases}$$

Introducing $\phi(r)$ as the smooth part of $\psi(r)$:  

$$\psi(r) = e^{-i\omega r_*} \phi(r),$$

the equation transforms to:

$$f(r)\phi''(r) + \left[f'(r) - 2i\omega\right]\phi'(r) - \frac{V(r)}{f(r)}\phi(r) = 0.$$  

To make coefficients of the equation finite for $r \in [0, +\infty)$, we define a new function:

$$F = \frac{\Delta}{r^4}$$

Then, we have:

$$F(r) \phi''(r) + \left( F'(r) + \frac{2 (r F - i \omega)}{r^2} \right) \phi'(r) - \left[ \frac{r F'(r)+2 F}{r^2} + \frac{l (l+1)}{r^4} \right] \phi = 0 $$

For numerical computation, compactify the radial coordinate using:

$$x \equiv 1 - \frac{r_+}{r}.$$

Here, $x \in [0, 1]$.
"""

# ╔═╡ bd2492fc-e012-42ea-bba5-3718ff10d6b4
md"""
**Input angular number $\ell$:**
"""

# ╔═╡ 9006db87-f4fb-42c5-b3ab-46e1434e3b42
@bind l NumberField(0:100, default=0)

# ╔═╡ 0ee959ae-25b1-4f4c-aa1d-1deb4b1204cc
md"""
**Input degree of the polynomial $n$:**
"""

# ╔═╡ 3a1ad375-c718-4465-95c6-d2f5bb9093df
@bind n NumberField(10:1000, default=40)

# ╔═╡ d9c28ff1-8b86-4027-9ff3-0557d9565f79
qnm = begin
	# domain x ∈ [0, 1]
	x_min = zero(TF)
	x_max = one(TF)
	dom = x_min .. x_max

	chebSpace = Chebyshev(dom)
	D = Derivative(chebSpace)
	conversion1to2 = Conversion(Ultraspherical(1, dom), Ultraspherical(2, dom))

	x = Fun(chebSpace)
	F = (r₊ * (r₊^3 / R^2 + r₊ * (x - 1)^2 + 2 * M * (x - 1)^3) + Q^2 * (x - 1)^4) / r₊^4
	dF = F'

	# generalized eigenvalue problem A ϕ = ω B ϕ

	# construct operators A
	c02 = (x - 1)^2 * F
	c01 = (x - 1)^2 * dF
	c00 = -2 * F + (x - 1) * (- l * (l + 1) * (x - 1) / r₊^2 + dF)
	A = c02 * D^2 + c01 * D + c00
	Am = @view(A[1:n, 1:n]) |> Matrix |> complex

	# construct matrix B
	c11 = 2im * (x - 1)^2 / r₊
	B = conversion1to2 * (c11 * D)
	Bm = @view(B[1:n, 1:n]) |> Matrix

	# enforce Dirichlet boundary conditions at x = 1
	Am[end, :] .= 1
	Bm[end, :] .= 0
	
	# Solve the generalized eigenvalue problem A ϕ = ω B ϕ
	qnm = eigvals!(Am, Bm; sortby=abs)
	filter(x -> !isnan(x) && imag(x) < 0, qnm)
end

# ╔═╡ abdb5714-4794-4435-ad92-0fd349fa7e77
Markdown.MD(Markdown.Admonition("warning", "TODO", [md"Ensure you perform a convergence test before publishing the results; this notebook is just a **tutorial**!"]))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ApproxFun = "28f2ccd6-bb30-5033-b560-165f7b14dc2f"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
ApproxFun = "~0.13.28"
LaTeXStrings = "~1.4.0"
PlutoUI = "~0.7.61"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "199c4e11031110f1214b6a9590fa9dca4ccd8f9e"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ApproxFun]]
deps = ["AbstractFFTs", "ApproxFunBase", "ApproxFunFourier", "ApproxFunOrthogonalPolynomials", "ApproxFunSingularities", "Calculus", "DomainSets", "FastTransforms", "LinearAlgebra", "RecipesBase", "Reexport", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "aaa2506645c1b7e1e7fdcf6182852de53d7a1bab"
uuid = "28f2ccd6-bb30-5033-b560-165f7b14dc2f"
version = "0.13.28"

    [deps.ApproxFun.extensions]
    ApproxFunDualNumbersExt = "DualNumbers"

    [deps.ApproxFun.weakdeps]
    DualNumbers = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"

[[deps.ApproxFunBase]]
deps = ["AbstractFFTs", "BandedMatrices", "BlockArrays", "BlockBandedMatrices", "Calculus", "Combinatorics", "DSP", "DomainSets", "FFTW", "FillArrays", "InfiniteArrays", "IntervalSets", "LazyArrays", "LinearAlgebra", "LowRankMatrices", "SparseArrays", "SpecialFunctions", "StaticArrays", "Statistics"]
git-tree-sha1 = "a007883ec340f24080d8fe4fc5518cabc54aea02"
uuid = "fbd15aa5-315a-5a7d-a8a4-24992e37be05"
version = "0.9.32"

    [deps.ApproxFunBase.extensions]
    ApproxFunBaseDualNumbersExt = "DualNumbers"
    ApproxFunBaseTestExt = "Test"

    [deps.ApproxFunBase.weakdeps]
    DualNumbers = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ApproxFunFourier]]
deps = ["AbstractFFTs", "ApproxFunBase", "BandedMatrices", "DomainSets", "FFTW", "FastTransforms", "InfiniteArrays", "IntervalSets", "LinearAlgebra", "Reexport", "StaticArrays"]
git-tree-sha1 = "0cb3a0f77266ae5e13f5ac155dc2082b7aea47c9"
uuid = "59844689-9c9d-51bf-9583-5b794ec66d30"
version = "0.3.31"

[[deps.ApproxFunOrthogonalPolynomials]]
deps = ["ApproxFunBase", "BandedMatrices", "BlockArrays", "BlockBandedMatrices", "DomainSets", "FastGaussQuadrature", "FastTransforms", "FillArrays", "HalfIntegers", "IntervalSets", "LinearAlgebra", "OddEvenIntegers", "Reexport", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "086fbeec1145f5b25d91ccbb8b7ce482b6ea3f3d"
uuid = "b70543e2-c0d9-56b8-a290-0d4d6d4de211"
version = "0.6.59"

    [deps.ApproxFunOrthogonalPolynomials.extensions]
    ApproxFunOrthogonalPolynomialsPolynomialsExt = "Polynomials"
    ApproxFunOrthogonalPolynomialsStaticExt = "Static"

    [deps.ApproxFunOrthogonalPolynomials.weakdeps]
    Polynomials = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
    Static = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"

[[deps.ApproxFunSingularities]]
deps = ["ApproxFunBase", "ApproxFunOrthogonalPolynomials", "BlockBandedMatrices", "DomainSets", "HalfIntegers", "IntervalSets", "LinearAlgebra", "OddEvenIntegers", "Reexport", "SpecialFunctions"]
git-tree-sha1 = "7061a7d8f22a40cca68216d9af91c03880ba91ac"
uuid = "f8fcb915-6b99-5be2-b79a-d6dbef8e6e7e"
version = "0.3.21"
weakdeps = ["StaticArrays"]

    [deps.ApproxFunSingularities.extensions]
    ApproxFunSingularitiesStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "4e25216b8fea1908a0ce0f5d87368587899f75be"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "bbc6688495b031d84610e227d46c35e17fdde5f5"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.9.1"
weakdeps = ["SparseArrays"]

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

[[deps.BlockArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "1ded9033f6067573314b27cd4b9ff01a1ba92cff"
uuid = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
version = "1.4.0"
weakdeps = ["BandedMatrices"]

    [deps.BlockArrays.extensions]
    BlockArraysBandedMatricesExt = "BandedMatrices"

[[deps.BlockBandedMatrices]]
deps = ["ArrayLayouts", "BandedMatrices", "BlockArrays", "FillArrays", "LinearAlgebra", "MatrixFactorizations"]
git-tree-sha1 = "4eef2d2793002ef8221fe561cc822eb252afa72f"
uuid = "ffab5731-97b5-5995-9138-79e8c1846df0"
version = "0.13.4"
weakdeps = ["SparseArrays"]

    [deps.BlockBandedMatrices.extensions]
    BlockBandedMatricesSparseArraysExt = "SparseArrays"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9cb23bbb1127eefb022b022481466c0f1127d430"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.2"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.DSP]]
deps = ["Bessels", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "60b9e054b777bb8bf0fae29b1c9d17f6876b2b44"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.8.2"

    [deps.DSP.extensions]
    OffsetArraysExt = "OffsetArrays"

    [deps.DSP.weakdeps]
    OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "a7e9f13f33652c533d49868a534bfb2050d1365f"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.15"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "7de7c78d681078f027389e067864a8d53bd7c3c9"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4d81ed14783ec49ce9f2e168208a12ce1815aa25"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+3"

[[deps.FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "fd923962364b645f3719855c88f7074413a6ad92"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "1.0.2"

[[deps.FastTransforms]]
deps = ["AbstractFFTs", "ArrayLayouts", "BandedMatrices", "FFTW", "FastGaussQuadrature", "FastTransforms_jll", "FillArrays", "GenericFFT", "LazyArrays", "Libdl", "LinearAlgebra", "RecurrenceRelationships", "Reexport", "SpecialFunctions", "ToeplitzMatrices"]
git-tree-sha1 = "e914dd1c91d1909f6584d61d767fd0c4f64fcba3"
uuid = "057dd010-8810-581a-b7be-e3fc3b93f78c"
version = "0.16.8"

[[deps.FastTransforms_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "FFTW_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl", "MPFR_jll", "OpenBLAS_jll"]
git-tree-sha1 = "efb41482692019ed03e0de67b9e48e88c0504e7d"
uuid = "34b6f7d7-08f9-5794-9e10-3819e4c7e49a"
version = "0.6.3+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

    [deps.FillArrays.weakdeps]
    PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GMP_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "781609d7-10c4-51f6-84f2-b8444358ff6d"
version = "6.3.0+0"

[[deps.GenericFFT]]
deps = ["AbstractFFTs", "FFTW", "LinearAlgebra", "Reexport"]
git-tree-sha1 = "1bc01f2ea9a0226a60723794ff86b8017739f5d9"
uuid = "a8297547-1b15-4a5a-a998-a2ac5f1cef28"
version = "0.1.6"

[[deps.HalfIntegers]]
git-tree-sha1 = "9c3149243abb5bc0bad0431d6c4fcac0f4443c7c"
uuid = "f0d1745a-41c9-11e9-1dd9-e5d34d218721"
version = "1.6.0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InfiniteArrays]]
deps = ["ArrayLayouts", "FillArrays", "Infinities", "LazyArrays", "LinearAlgebra"]
git-tree-sha1 = "7e3e697ea3d74435d7d396c75d1356f3bcc6e9cf"
uuid = "4858937d-0d70-526a-a4dd-2d5cb5dd786c"
version = "0.15.5"
weakdeps = ["BandedMatrices", "BlockArrays", "BlockBandedMatrices", "DSP", "Statistics"]

    [deps.InfiniteArrays.extensions]
    InfiniteArraysBandedMatricesExt = "BandedMatrices"
    InfiniteArraysBlockArraysExt = "BlockArrays"
    InfiniteArraysBlockBandedMatricesExt = "BlockBandedMatrices"
    InfiniteArraysDSPExt = "DSP"
    InfiniteArraysStatisticsExt = "Statistics"

[[deps.Infinities]]
git-tree-sha1 = "437cd9b81b649574582507d8b9ca3ffc981f2fc2"
uuid = "e1ba4f0e-776d-440f-acd9-e1d2e9742647"
version = "0.1.11"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "0f14a5456bdc6b9731a5682f439a672750a09e48"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.0.4+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "8c739d0eee1d00f3411dcc8c28ba732f84480635"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.6.0"
weakdeps = ["BandedMatrices", "BlockArrays", "BlockBandedMatrices", "StaticArrays"]

    [deps.LazyArrays.extensions]
    LazyArraysBandedMatricesExt = "BandedMatrices"
    LazyArraysBlockArraysExt = "BlockArrays"
    LazyArraysBlockBandedMatricesExt = "BlockBandedMatrices"
    LazyArraysStaticArraysExt = "StaticArrays"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LowRankMatrices]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7c8664b2f3d5c3d9b77605c03d53b18813e79b0f"
uuid = "e65ccdef-c354-471a-8090-89bec1c20ec3"
version = "1.0.1"
weakdeps = ["FillArrays"]

    [deps.LowRankMatrices.extensions]
    LowRankMatricesFillArraysExt = "FillArrays"

[[deps.MIMEs]]
git-tree-sha1 = "1833212fd6f580c20d4291da9c1b4e8a655b128e"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.0.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "5de60bc6cb3899cd318d80d627560fae2e2d99ae"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.0.1+1"

[[deps.MPFR_jll]]
deps = ["Artifacts", "GMP_jll", "Libdl"]
uuid = "3a97d323-0669-5f0c-9066-3539efd106a3"
version = "4.2.1+0"

[[deps.MacroTools]]
git-tree-sha1 = "72aebe0b5051e5143a079a4685a46da330a40472"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.15"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "16a726dba99685d9e94c8d0a8f655383121fc608"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "3.0.1"
weakdeps = ["BandedMatrices"]

    [deps.MatrixFactorizations.extensions]
    MatrixFactorizationsBandedMatricesExt = "BandedMatrices"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OddEvenIntegers]]
git-tree-sha1 = "9980c824d61373cd473149480d3ef59119f382bf"
uuid = "8d37c425-f37a-4ca2-9b9d-a61bc06559d2"
version = "0.1.12"
weakdeps = ["HalfIntegers"]

    [deps.OddEvenIntegers.extensions]
    OddEvenIntegersHalfIntegersExt = "HalfIntegers"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "7e71a55b87222942f0f9337be62e26b1f103d3e4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.61"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "0973615c3239b1b0d173b77befdada6deb5aa9d8"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.17"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecurrenceRelationships]]
git-tree-sha1 = "1a268b6c7f2eacf86dc345d85a88d7120c4f8103"
uuid = "807425ed-42ea-44d6-a357-6771516d7b2c"
version = "0.1.1"
weakdeps = ["FillArrays", "LazyArrays", "LinearAlgebra"]

    [deps.RecurrenceRelationships.extensions]
    RecurrenceRelationshipsFillArraysExt = "FillArrays"
    RecurrenceRelationshipsLazyArraysExt = "LazyArrays"
    RecurrenceRelationshipsLinearAlgebraExt = "LinearAlgebra"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "64cca0c26b4f31ba18f13f6c12af7c85f478cfde"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "e3be13f448a43610f978d29b7adf78c76022467a"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.12"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.ToeplitzMatrices]]
deps = ["AbstractFFTs", "DSP", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "338d725bd62115be4ba7ffa891d85654e0bfb1a1"
uuid = "c751599d-da0a-543b-9d20-d0a503d91d24"
version = "0.8.5"

    [deps.ToeplitzMatrices.extensions]
    ToeplitzMatricesStatsBaseExt = "StatsBase"

    [deps.ToeplitzMatrices.weakdeps]
    StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═fdf0edd0-efbc-11ef-1dea-e32569d65705
# ╠═64dcba7b-daea-488c-bc89-e729cc7537e7
# ╟─e57eb498-14bc-4606-89d0-3ddbaa20c5b5
# ╠═6327b564-6066-45b4-b188-095e3761b96d
# ╟─b1db4296-da97-4abe-9303-b77878f212ac
# ╟─8d2ee7b2-e5e9-402a-95b0-b5c53db84c67
# ╠═4d0cd86e-8cb4-4c40-a202-6513b6257892
# ╟─24060d90-60f9-4dd6-8206-a28514fc2ab1
# ╠═63dc53d0-0c1c-432a-9a0e-07502660eb6d
# ╠═a5869853-874c-4779-866f-ae5e7af251a8
# ╟─053e3de3-8ea4-409c-aa00-5d2dc07d8f14
# ╠═07200c40-621c-46f0-ae95-58a3b4755648
# ╟─33b34c18-883c-46e2-bc84-80f5e411d382
# ╠═e687ba54-efd4-4d0d-967c-eb20160a8c51
# ╠═b766f294-5c18-40ce-8f92-696c1dcaaacb
# ╟─7f6a4b13-c252-47e9-b01c-5db1a7710985
# ╠═a5dc2c36-0c26-4cfe-a03a-137aabb47467
# ╟─9a09cfdf-98ae-4c73-943a-91472c0798e8
# ╟─bd2492fc-e012-42ea-bba5-3718ff10d6b4
# ╠═9006db87-f4fb-42c5-b3ab-46e1434e3b42
# ╟─0ee959ae-25b1-4f4c-aa1d-1deb4b1204cc
# ╠═3a1ad375-c718-4465-95c6-d2f5bb9093df
# ╠═d9c28ff1-8b86-4027-9ff3-0557d9565f79
# ╠═abdb5714-4794-4435-ad92-0fd349fa7e77
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
