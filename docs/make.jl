using Documenter, PortHamiltonianBenchmarkSystems

push!(LOAD_PATH, "../src/")

makedocs(
    sitename = "PortHamiltonianBenchmarkSystems",
    pages = [
        "Home" => "index.md",
        # "Contribution" => "contribution.md",
        "Benchmark Systems" =>
            ["GugercinMSDChain.md", "PoroModel.md", "RclCircuits.md", "DampedWaveNet.md"],
    ],
)

deploydocs(
    repo = "github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl.git",
    versions = nothing,
)
