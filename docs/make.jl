using GridFit
using Documenter

DocMeta.setdocmeta!(GridFit, :DocTestSetup, :(using GridFit); recursive=true)

makedocs(;
    modules=[GridFit],
    authors="azraq27 <bill.gross@me.com> and contributors",
    repo="https://github.com/azraq27/GridFit.jl/blob/{commit}{path}#{line}",
    sitename="GridFit.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://azraq27.github.io/GridFit.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/azraq27/GridFit.jl",
    devbranch="main",
)
