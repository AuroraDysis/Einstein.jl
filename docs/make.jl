using PDESuite
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(PDESuite, :DocTestSetup, :(using PDESuite); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "src/refs.bib"); style=:alpha)

makedocs(;
    modules=[PDESuite],
    authors="Zhen Zhong <auroradysis@gmail.com> and contributors",
    sitename="PDESuite.jl",
    format=Documenter.HTML(;
        canonical="https://auroradysis.github.io/PDESuite.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Chebyshev Suite" => "cheb.md",
        "Finite Difference Suite" => "fdm.md",
    ],
    plugins=[bib],
)

deploydocs(; repo="github.com/AuroraDysis/PDESuite.jl", devbranch="main")
