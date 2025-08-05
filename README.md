# PortHamiltonianBenchmarkSystems.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://algopaul.github.io/PortHamiltonianBenchmarkSystems.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://algopaul.github.io/PortHamiltonianBenchmarkSystems.jl/dev/)
[![Build Status](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/actions/workflows/CI.yml)
[![codecov](http://codecov.io/github/Algopaul/PortHamiltonianBenchmarkSystems.jl/coverage.svg?branch=main)](http://codecov.io/github/Algopaul/PortHamiltonianBenchmarkSystems.jl?branch=main)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

[PortHamiltonianBenchmarkSystems](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/) is a collection of port-Hamiltonian systems, that can be used as benchmarks for simulation, control, and model-order reduction algorithms. We feature constructors for a wide range of linear, nonlinear, ODE, and DAE systems, as well as several default parameter sets for each. Detailed descriptions and derivations of the featured systems can be found in the [documentation](https://algopaul.github.io/PortHamiltonianBenchmarkSystems.jl/).

If you want to take advantage of this benchmark collection programming languages other than julia, you can:
- Generate any desired system in Julia and save the matrices in a format of your choosing (see [JuliaIO](https://github.com/JuliaIO)),
- Generate MAT-files for any default parameter set, using our [Command-Line Interface](https://github.com/Algopaul/PortHamiltonianBenchmarkSystemsCLI.jl).

## Installation and Usage

Install with the Julia package manager [Pkg](https://pkgdocs.julialang.org/):
```julia
pkg> add PortHamiltonianBenchmarkSystems # Press ']' to enter the Pkg REPL mode.
```
or in the [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/):
```julia
julia> using Pkg; Pkg.add("PortHamiltonianBenchmarkSystems")
```
To generate one of the systems, e.g. a mass-spring-damper chain with the parameters from [Gugercin2012](https://doi.org/10.1016/j.automatica.2012.05.052), type:
```julia
using PortHamiltonianBenchmarkSystems
config = SingleMSDConfig("Gugercin")
J, R, Q, B = construct_system(config)
```
Naturally, we may also specify the parameters ourselves:
```julia
config = SingleMSDConfig(n_cells = 10, io_dim = 2, c = 1.0, m = 4.0, k = 5.0)
```
If you need the system matrices in standard port-Hamiltonian form, type:
```julia
system = PHSystem(config)
(; E, J, R, Q, G, P, S, N) = system
```
Docstrings for the types and methods shown above can be accessed in the Julia REPL by typing `?` and then name of the type or method.

## How to Contribute

This benchmark collection is driven by the active support of the port-Hamiltonian community. If your research has lead to port-Hamiltonian models that may be relevant for this collection, we would be happy to include them. If you wish to contribute to the project directly, please consult our [Contribution](@ref) page, fork our [Git repository](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/) and issue a pull request. Alternatively, feel free to contact us via [e-mail](mailto:paul.schwerdtner@nyu.edu) to discuss our potential collaboration. We are happy to receive reference implementations in other languages and reimplement them in Julia for this package.
