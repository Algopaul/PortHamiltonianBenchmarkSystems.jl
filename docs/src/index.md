# PortHamiltonianBenchmarkSystems

This is a collection of port-Hamiltonian benchmark systems that can be used to test algorithms designed for simulation, control, or model-order reduction of port-Hamiltonian systems. We feature a wide spectrum of linear, nonlinear, ODE and DAE systems.

**PortHamiltonianBenchmarkSystems** is a `julia`-Package that also supports `MATLAB` or `Python` users by providing download links to all benchmark examples. These can be found in the [WIKI](https://github.com/Amanibus/PortHamiltonianBenchmarkSystems/wiki). Data is created via a CI pipeline, so users from ``julia``, ``MATLAB``, and ``Python`` can be sure to work with the same systems.

## Installation and Usage

In a ``julia`` instance, execute
```julia
using Pkg
Pkg.add(url="https://github.com/Algopaul/PortHamiltonianBenchmarkSystems/")
```
to install **PortHamiltonianBenchmarkSystems** and gain access to all BenchmarkExamples directly within julia. To load the one of the benchmark systems simply type
```julia
using PortHamiltonianBenchmarkSystems
J, R, Q, B = gugercin_pH_msd_chain()
```
to load the model used in [Gugercin2012](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems/blob/7c7e588f9bd67ba4a5c67ac37768c9c43021e6e6/bibliography.tex#L9-L17) to test the effectiveness of PH-IRKA. To learn about tunable options for the example, simply execute
```julia
?gugercin_pH_msd_chain
```
in the ``julia`` REPL.

## Contribution

If you want to contribute your PH-model to the collection, ...

