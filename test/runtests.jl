using PortHamiltonianBenchmarkSystems
using LinearAlgebra
using Test

@testset "PortHamiltonianBenchmarkSystems" verbose = true begin

    include("./SingleMSDChainTests.jl")
    include("./PoroElasticityModelTests.jl")
    @testset "RCL-Ladders" begin end

end
