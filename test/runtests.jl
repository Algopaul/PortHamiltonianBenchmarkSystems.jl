using PortHamiltonianBenchmarkSystems
using LinearAlgebra
using Test

@testset "PortHamiltonianBenchmarkSystems" verbose = true begin
    include("./SingleMSDChainTests.jl")
    include("./PoroElasticityModelTests.jl")
    include("./RCLLaddersTests.jl")
    include("./DampedWaveNetTests.jl")
end
