using PortHamiltonianBenchmarkSystems
using LinearAlgebra
using SparseArrays
using Test, Aqua

@testset "PortHamiltonianBenchmarkSystems" verbose = true begin
    Aqua.test_all(PortHamiltonianBenchmarkSystems)

    include("./SingleMSDChainTests.jl")
    include("./PoroElasticityModelTests.jl")
    include("./RCLLadderODETests.jl")
    include("./RCLLaddersTests.jl")
    include("./DampedWaveNetTests.jl")
end
