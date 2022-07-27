# PortHamiltonianBenchmarkSystems

## About

**PortHamiltonianBenchmarkSystems.jl** is a collection of port-Hamiltonian system constructors, that can be used as benchmarks for simulation, control and model-order reduction algorithms. We feature a wide spectrum of linear, nonlinear, ODE and DAE systems, all of which can be generated with full choice of parameters. If you do not plan to use Julia as a main language for your project, but still want to take advantage of this package, you can download `mat`-files for each system in our collection at link. Alternatively, you could generate your desired system in Julia and save the matrices as `mat`-files, using [MAT.jl](https://github.com/JuliaIO/MAT.jl).

## Installation and Usage

To install **PortHamiltonianBenchmarkSystems** and gain access to all benchmark systems, run the following in the Julia REPL:
```julia
using Pkg
Pkg.add(url="https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/")
```
To generate one of the systems, e.g. a single mass-spring-damper chain with the parameters used in [Gugercin2012](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems/blob/7c7e588f9bd67ba4a5c67ac37768c9c43021e6e6/bibliography.tex#L9-L17), simply type:
```julia
using PortHamiltonianBenchmarkSystems

config = MSDChainConfig("Gugercin")
J,R,Q,B = construct_system(config)
```
Naturally, we may change any parameters we wish, or specify all of them ourselves:
```julia
config = MSDChainConfig("Gugercin")
config.n = 20
config.k = 5.

n,d,c,m,k = (20,2,1.,4.,5.)
config = MSDChainConfig(n,d,c,m,k)
```
If you need the system matrices in standard port-Hamiltonian form, simply type:
```julia
system = PHSystem(config)
```
The full documentation of the benchmark system can be loaded by typing:
```julia
?SingleMSDChain
```

## Contribution

This benchmark collection is driven by the active support of the port-Hamiltonian community. If your research has lead to port-Hamiltonian models that may be relevant for this collection, we would be happy to include them. All entries consist of code for generating system matrices, a comprehensive set of tests and a documentation page, showing the derivation of the system and a code reference. Exact instructions for each of these parts are provided in the sections below. To add your benchmark model to the collection, either fork our repository, add your entry and issue a pull request, or send us your contribution directly via [E-Mail](mailto:schwerdt@math.tu-berlin.de).
