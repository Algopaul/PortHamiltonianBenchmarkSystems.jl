using Documenter, PortHamiltonianBenchmarkSystems

push!(LOAD_PATH, "../src/")

makedocs(
    sitename = "PortHamiltonianBenchmarkSystems",
    pages = [
        "Home" => "index.md",
        "Contribution" => "Contribution.md",
        "Benchmark Systems" =>
            ["SingleMSDChain.md", "PoroModel.md", "RclCircuits.md", "DampedWaveNet.md", "Elasticity2DAFW.md"],
    ],
)

deploydocs(
    repo = "github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl.git",
    versions = nothing,
)
