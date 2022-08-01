# PortHamiltonianBenchmarkSystems

## About

[PortHamiltonianBenchmarkSystems](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/) is a collection of port-Hamiltonian systems, that can be used as benchmarks for simulation, control, and model-order reduction algorithms. We feature constructors for a wide range of linear, nonlinear, ODE, and DAE systems, as well as several default parameter sets for each. Detailed descriptions of each featured system can be found in the sidebar, under 'Benchmark Systems'.

!!! note
    If you do not intend to use Julia as a main language in your research, but still want to take advantage of this package, you can download MAT-files for all default configurations of each system from our [Zenodo](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/). Alternatively, you could generate your desired systems in Julia and save them in a format of your choosing (e.g. as MAT-files, using [MAT.jl](https://github.com/JuliaIO/MAT.jl)). Additionally, we provide a command-line-interface for the benchmark collection [here](https://github.com/Algopaul/PortHamiltonianBenchmarkSystemsCLI.jl).

## Installation and Usage

To install [PortHamiltonianBenchmarkSystems](https://algopaul.github.io/PortHamiltonianBenchmarkSystems.jl/), run the following commands in the Julia REPL:
```julia
using Pkg
Pkg.add(url="https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/")
```
To generate one of the systems, e.g. a mass-spring-damper chain with the parameters from [Gugercin2012](), type:
```julia
using PortHamiltonianBenchmarkSystems

config = SingleMSDConfig("Gugercin")
J, R, Q, B = construct_system(config)
```
Naturally, we may also specify the parameters ourselves.
```julia
config = SingleMSDConfig(10, 2, 1.0, 4.0, 5.0)
```

If you need the system matrices in standard port-Hamiltonian form, type:
```julia
system = PHSystem(config)
E, J, R, Q, G, P, S, N = @unpack system
```
The docstrings for the constructors and methods shown above can be pulled up in the Julia REPL as follows. Documentation is accessed in the julia REPL by typing `?` followed by the name of the function or struct.

For a detailed description of the system in question, an overview of its discretization and a code reference, consult the [Single MSD Chain](@ref) page.

## How to Contribute

This benchmark collection is driven by the active support of the port-Hamiltonian community. If your research has lead to port-Hamiltonian models that may be relevant for this collection, we would be happy to include them. If you wish to contribute to the project directly, please consult our [Contribution](@ref) page, fork our [git repository](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/) and issue a pull request. Alternatively, feel free to contact us via [e-mail](mailto:schwerdt@math.tu-berlin.de) to discuss our potential collaboration. We are happy to receive reference implementations in other languages and reimplement them in Julia for this package.
