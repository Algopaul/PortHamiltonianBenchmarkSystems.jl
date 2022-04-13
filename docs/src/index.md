# PortHamiltonianBenchmarkSystems

## About

This is a collection of port-Hamiltonian benchmark systems that can be used to test algorithms designed for simulation, control, or model-order reduction of port-Hamiltonian systems. We feature a wide spectrum of linear, nonlinear, ODE and DAE systems. Benchmark systems are loaded via functions that also allow to set free modeling parameters (such as the stiffness in a mass spring damper chain).

**PortHamiltonianBenchmarkSystems** is a julia-Package that also supports MATLAB or python users by providing download links to all benchmark examples.

## Installation and Usage

In the julia REPL, execute
```julia
using Pkg
Pkg.add(url="https://github.com/Algopaul/PortHamiltonianBenchmarkSystems/")
```
to install **PortHamiltonianBenchmarkSystems** and gain access to all
BenchmarkExamples directly within julia. To load the one of the benchmark
systems simply type
```julia
using PortHamiltonianBenchmarkSystems
J, R, Q, B = gugercin_pH_msd_chain()
```
to load the model used in [Gugercin2012](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems/blob/7c7e588f9bd67ba4a5c67ac37768c9c43021e6e6/bibliography.tex#L9-L17). To learn about tunable options for this benchmark example, simply load the documentation using the following call.
```julia
?gugercin_pH_msd_chain
```

!!! note

    If you do not plan to use julia as main language for your project but still
    want to take advantage of this package, you can download `mat`-files for
    each benchmark example in our collection. The links can be found in this
    documentation at the corresponding example page. However, in this way you
    cannot configure the parameters yourself. For that, you can first generate
    the system matrices in julia and then use [MAT](https://github.com/JuliaIO/MAT.jl) package to store them as
    `mat`-file and use that in your preferred programming environment.

## Contribution

This benchmark collection is driven by the active support of the port-Hamiltonian community. If your research has lead to port-Hamiltonian models that may be relevant for this collection, we are happy to include it in this collection. Each benchmark system consists of code that generates the corresponding system matrices as well as a documentation that explains the origin, potential application, and special features of the model. For reference, check out our first example pages for a [port-Hamiltonian mass spring damper chain](./GugercinMSDChain.md) and the corresponding [code](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems/blob/58e925c50836958a83141ae987b0b5ace4be953c/src/PortHamiltonianBenchmarkSystems.jl#L25).

To add your example, you can either simply fork our repository, add your example generation code and documentation, and issue a pull request or send us your code and a markdown file (using [``\KaTeX``](https://katex.org/) for math expressions) to [Paul Schwerdtner](mailto:schwerdt@math.tu-berlin.de).
!!! warning

    If your benchmark example requires large files for constructions, do not try to add them to this git repository. Instead we recommend to use Zenodo for storing the large files and downloading them on request. An example for this can be found [here](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems/blob/58e925c50836958a83141ae987b0b5ace4be953c/src/PortHamiltonianBenchmarkSystems.jl#L107).
