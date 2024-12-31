using GRSuite
using Documenter

DocMeta.setdocmeta!(GRSuite, :DocTestSetup, :(using GRSuite); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "src/refs.bib"))

makedocs(;
    modules=[GRSuite],
    authors="Zhen Zhong <auroradysis@gmail.com> and contributors",
    sitename="GRSuite.jl",
    format=Documenter.HTML(;
        canonical="https://auroradysis.github.io/GRSuite.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
    plugins=[bib],
)

deploydocs(; repo="github.com/AuroraDysis/GRSuite.jl", devbranch="main")
