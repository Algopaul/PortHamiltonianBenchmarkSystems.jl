# Contribution

## Modus Operandi
This benchmark collection is driven by the active support of the port-Hamiltonian community. If your research has lead to port-Hamiltonian models that may be relevant for this collection, we would be happy to include them. If you wish to contribute to the project directly, please consult this page, fork our [git repository](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/) and issue a pull request. Alternatively, feel free to contact us via [e-mail](mailto:schwerdt@math.tu-berlin.de) to discuss our potential collaboration. We are happy to receive reference implementations in other languages and reimplement them in Julia for this package.

All entries in the collection consist of code for generating the system matrices, a comprehensive set of tests and a documentation page. Exact instructions for each of these parts are provided in the sections below. Some of the design choices may still be subject to change, but efforts are currently being made to settle on a final structure. Any subsequently required changes to community contributions will be handled by us.

## Code
This package consists of a single module, containing the following six elements:
- `PHSystem`: parametric composite type for storing system matrices in port-Hamiltonian form, with a single input validating internal constructor;
- `BenchmarkConfig`: abstract type for grouping all `<Model>Config` types;
- `<Model>Config <: BenchmarkConfig`: (parametric) composite type containing all parameters for some model, with a single input validating internal constructor;
- `<Model>Config`: external constructor providing several default `<Model>Config` instances based on some identifier
- `construct_system`: method for constructing system matrices in 'natural' form based on some `<Model>Config` instance;
- `PHSystem`: external constructor for constructing system matrices in port-Hamltonian form.
The last four elements are repeated for each benchmark model and are stored together in `src/<Model>.jl`. All `<Model>.jl` files are included in `src/PortHamiltonianBenchmarkSystems.jl`, which also contains the type declarations for `PHSystem` and `BenchmarkConfig`.

To contribute to the code, simply add a file based on the example below at `/src/<Model>.jl` and add a corresponding `include` statement to `src/PortHamiltonianBenchmarkSystems.jl`. As shown below, the docstrings should be written in `Markdown` format and any default `<Model>Config` parameters should be stored on our [Zenodo](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/) and retreived as Julia artifacts when needed.
```julia
export RandLinConfig

"""
Composite type describing a linear port-Hamiltonian system, where all independent matrix entries are randomly chosen between 0 and 1.
# Arguments
- n_x: number of state variables
- n_p: number of input and output ports
"""
struct RandLinConfig <: BenchmarkConfig
    n_x::Int64
    n_p::Int64

    RandLinConfig(n_x::Int64, n_p::Int64)
        #Validate parameters
        @assert (n_x, n_p) .> 0 "Number of state variables and ports must be larger than 0"
        @assert n_p  <= n_x "???"

        return new(n_x, n_p)
    end
end

"""
External constructor for retrieving default RandLinConfig instances.
# Arguments
- `id`: identifier for default parameter set
- `n_x`: override for ``n_x``
- `n_p`: override for ``n_p``
"""
function RandLinConfig(id::String; n_x = nothing, n_p = nothing)
    params = matread(artifact"pH_RandLinConfig_" * id)    
    n_x == nothing ? n_x = params["n_x"] :
    n_p == nothing ? n_p = params["n_p"] :

    return RandLinConfig(n_x, n_p)
end

"""
Method for constructing the system matries in natural form.
# Arguments
- `config`: `RandLinConfig` instance
# Output
- `system`: Named tuple containing the system matrices in 'natural' form
"""
function construct_system(config::RandLinConfig)
    M1 = rand(config.n_x, config.n_x)
    M2 = rand(config.n_x, config.n_x)
    M3 = rand(config.n_x, config.n_p)

    J = (M1 - M1')/sqrt(2)
    R = (M2 * M2')/sqrt(config.n_x)
    G = (M3 - M3')/sqrt(3)

    return (J = J, R = R, G = G)
end

function PHSystem(config:RandLinConfig)
    J, R, G = construct_system(config)

    E = sparse(1.0I,size(J)...)
    Q = sparse(1.0I,size(J)...)
    P = spzeros(size(G)...)
    N = spzeros(config.n_p, config.n_p)
    S = spzeros(size(N)...)

    return PHSystem(E, J, R, Q, G, P, S, N)
end
```

## Tests
To guarantee that all merged code is in a working state, automatic test pipelines have been set up for this repository. The test are run through GitHub Actions, using the `Test.jl` package. 

Every benchmark model has its own `/test/<Model>Tests.jl` file, which is included in `/test/runtests.jl`. While the goal is to achieve near total code coverage, it is not necessary to 'test everything'. We suggest to at least test the following for a variety of configurations:
- System matrix sizes;
- System transfer functions;
- Correspondence between 'natural' and port-Hamiltonian system matrices.

To contribute to the tests, add a file based on the template below at `/test/<Model>Tests.jl` and add a corresponding `include` statement to `/test/runtests.jl`. Please refer to the [current tests](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/tree/main/test) for more detailed examples.
```julia
@testset "RandLin" begin
    @test #Test 1
    @test #Test ...
    @test #Test N
end
```

## Documentation
The documentation for this package is built using `Documenter.jl`. The `/docs/make.jl` script uses the Markdown files in `/docs/src` and the images in `/docs/src/assets`to build a documentation webpage in `/docs/build`. Equations are rendered using [``\KaTeX``](https://katex.org/), which is invoked in Markdown by the `math` environment. The webpage can be loaded locally by running `make.jl` and then `LiveServer.serve(dir=/docs/build`). The documentation should then be accessible from the returned `http://localhost` port.

Each benchmark model is documented in a separate file, containing the following sections:
- `Description`: mathematical description of the model;
- `Discretization`: (if applicable) detailed description of the discretization procedure and conversion to port-Hamiltonian form;
- `Interface`: section importing the docstrings from the corresponding `<Model>.jl` file;
- `References`: reference section in `BibTeX` format.

To contribute to the documentation, add a file based on the example below at `/docs/src/<Model>.md` and add `"<Model>.md"` to the `"Benchmark Systems"` list in `/docs/make.jl`. Please use the current 'Benchmark System' pages for stylistic reference.
````markdown
# Random Linear pH-System

## Description
This benchmark is a linear port-Hamiltonian system:
```math
\begin{align*}
    E\dot{x} &= (J-R)Qx + (G-P)u\\
    y &= (G+P)^HQx + (S+N)u
\end{align*}
```
where ``E,Q=I``, ``P,S,N=0`` and ``J,\ R,\ G`` are random dense matrices of the correct structure (R positive semi-definite, J, G skew symmetric) with mean 0 and variance 1.

## Discretization
The system is discrete a priori.

## Interface
```@docs
RandLinConfig
```
```@docs
RandLinConfig(id::String)
```
```@docs
construct_system(config::RandLinConfig)
```

## References
```LaTeX

```
````
