using Documenter, PortHamiltonianBenchmarkSystems

push!(LOAD_PATH, "../src/")
makedocs(sitename="PortHamiltonianBenchmarkSystems Documentation")
deploydocs(
  repo = "github.com/Algopaul/PortHamiltonianBenchmarkSystems.git",
  versions = nothing
)
