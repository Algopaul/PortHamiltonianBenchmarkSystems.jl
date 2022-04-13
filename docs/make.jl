using Documenter, PortHamiltonianBenchmarkSystems

push!(LOAD_PATH, "../src/")
makedocs(
 sitename="PortHamiltonianBenchmarkSystems",
 pages = [
          "Home" => "index.md",
          "Benchmark Systems" => [
                                  "GugercinMSDChain.md", "PoroModel.md"
            ]
         ]
)
deploydocs(
  repo = "github.com/Algopaul/PortHamiltonianBenchmarkSystems.git",
  versions = nothing
)
