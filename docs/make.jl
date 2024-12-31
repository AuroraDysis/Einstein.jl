using GRSuite
using Documenter

DocMeta.setdocmeta!(GRSuite, :DocTestSetup, :(using GRSuite); recursive=true)

makedocs(;
    modules=[GRSuite],
    authors="Zhen Zhong <auroradysis@gmail.com> and contributors",
    sitename="GRSuite.jl",
    format=Documenter.HTML(;
        canonical="https://AuroraDysis.github.io/GRSuite.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/AuroraDysis/GRSuite.jl",
    devbranch="main",
)
