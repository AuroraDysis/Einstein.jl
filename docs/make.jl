using Einstein
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(Einstein, :DocTestSetup, :(using Einstein); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "src/refs.bib"); style=:alpha)

makedocs(;
    modules=[Einstein],
    authors="Zhen Zhong <auroradysis@gmail.com> and contributors",
    sitename="Einstein.jl",
    format=Documenter.HTML(;
        canonical="https://auroradysis.github.io/Einstein.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Chebyshev Suite" => "cheb.md",
        "Finite Difference Suite" => "fdm.md",
        "QNM Suite" => "qnm.md",
        "Utilities" => "utils.md",
    ],
    plugins=[bib],
)

deploydocs(; repo="github.com/AuroraDysis/Einstein.jl", devbranch="main")
