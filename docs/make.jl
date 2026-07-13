using CorePeriphery
using Documenter

makedocs(;
    modules=[CorePeriphery],
    checkdocs=:exports,
    authors="Simone Santoni",
    repo="https://github.com/simoneSantoni/CorePeriphery.jl/blob/{commit}{path}#{line}",
    sitename="CorePeriphery.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://simoneSantoni.github.io/CorePeriphery.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Algorithms" => "algorithms.md",
        "Contracts" => "contracts.md",
        "Performance" => "performance.md",
        "Development Record" => "development.md",
        "Migration Guide" => "migration.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/simoneSantoni/CorePeriphery.jl",
    devbranch="main",
)
