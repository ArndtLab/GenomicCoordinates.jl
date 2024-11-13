using GenomicCoordinates
using Documenter

DocMeta.setdocmeta!(GenomicCoordinates, :DocTestSetup, :(using GenomicCoordinates); recursive=true)

makedocs(;
    modules=[GenomicCoordinates],
    authors="Peter Arndt <arndt@molgen.mpg.de> and contributors",
    sitename="GenomicCoordinates.jl",
    format=Documenter.HTML(;
        canonical="https://ArndtLab.github.io/GenomicCoordinates.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    warnonly=[:missing_docs],
)

deploydocs(;
    repo="github.com/ArndtLab/GenomicCoordinates.jl",
    devbranch="main",
)
