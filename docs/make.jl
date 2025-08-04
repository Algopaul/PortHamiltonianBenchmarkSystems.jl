using PortHamiltonianBenchmarkSystems
using Documenter, DocumenterCitations

DocMeta.setdocmeta!(PortHamiltonianBenchmarkSystems, :DocTestSetup, :(using PortHamiltonianBenchmarkSystems); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "..", "CITATION.bib"))

makedocs(;
    modules=[PortHamiltonianBenchmarkSystems],
    sitename="PortHamiltonianBenchmarkSystems.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Algopaul.github.io/PortHamiltonianBenchmarkSystems.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Contribution" => "Contribution.md",
        "Benchmark Systems" =>
            ["SingleMSDChain.md", "PoroModel.md", "RCLLadderODE.md", "RclCircuits.md", "DampedWaveNet.md", "Elasticity2DAFW.md", "HeatModel.md", "LosslessWave.md"],
        "References" => "References.md",
    ],
    plugins=[bib],
)

deploydocs(
    repo = "github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl.git",
    devbranch = "main"
)
