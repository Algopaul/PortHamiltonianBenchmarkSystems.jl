# PortHamiltonianBenchmarkSystems.jl

## About

[PortHamiltonianBenchmarkSystems.jl](https://algopaul.github.io/PortHamiltonianBenchmarkSystems.jl/) is a collection of port-Hamiltonian system constructors, that can be used as benchmarks for simulation, control and model-order reduction algorithms. We feature a wide spectrum of linear, nonlinear, ODE and DAE systems, with full choice of parameters, as well as many default configurations.

!!! note
	If you do not intend to use Julia as a main language in your research or project, but still want to take advantage of this package, you can download MAT-files for various configurations of each system from our [Zenodo](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/). Alternatively, you could generate your desired systems in Julia and save them in a format of your choosing (e.g. as MAT-files, using [MAT.jl](https://github.com/JuliaIO/MAT.jl)).

## Installation and Usage

To install [PortHamiltonianBenchmarkSystems.jl](https://algopaul.github.io/PortHamiltonianBenchmarkSystems.jl/), run the following commands in the Julia REPL:
```julia
(@v1.7) pkg> add "https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/"
```
To generate one of the systems, e.g. a mass-spring-damper chain with the parameters from [Gugercin2012](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems/blob/7c7e588f9bd67ba4a5c67ac37768c9c43021e6e6/bibliography.tex#L9-L17), type:
```julia
using PortHamiltonianBenchmarkSystems

config = SingleMSDConfig("Gugercin")
J, R, Q, B = construct_system(config)
```
Naturally, we may modify any parameters we wish, or specify all of them ourselves:
```julia
config = SingleMSDConfig("Gugercin", k = 5.0)
config = SingleMSDConfig(10, 2, 1.0, 4.0, 5.0)
```
If you need the system matrices in standard port-Hamiltonian form, type:
```julia
system = PHSystem(config)
E, J, R, Q, G, P, S, N = @unpack system
```
The docstrings for the constructors and methods shown above can be pulled up in the Julia REPL as follows:
```julia
help?> SingleMSDChain
help?> construct_system(config::SingleMSDConfig)
help?> PHSystem(config::SingleMSDConfig)
```
For a detailed description of the system in question, an overview of its discretization and a code reference, consult the [Single MSD Chain](@ref) page.

## How to Contribute

This benchmark collection is driven by the active support of the port-Hamiltonian community. If your research has lead to port-Hamiltonian models that may be relevant for this collection, we would be happy to include them. If you wish to contribute to the project directly, please consult our [Contribution](@ref) page, fork our [git repository](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/) and issue a pull request. Alternatively, feel free to contact us via [e-mail](mailto:schwerdt@math.tu-berlin.de) to discuss our potential collaboration. We are happy to receive reference implementations in other languages and reimplement them in Julia for this package.
