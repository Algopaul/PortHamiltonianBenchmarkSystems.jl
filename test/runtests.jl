using PortHamiltonianBenchmarkSystems
using LinearAlgebra
using Test

@testset "PortHamiltonianBenchmarkSystems" verbose = true begin

    include("./SingleMSDChainTests.jl")

    @testset "PoroElasticityModel" begin end

    @testset "RCL-Ladders" begin end

end
