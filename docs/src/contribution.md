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
using Random

export RandLinConfig

"""
Composite type describing a linear port-Hamiltonian system, where all independent matrix entries are randomly chosen between 0 and 1.
# Arguments
- n_x: number of state variables
- n_u: number of inputs
"""
struct RandLinConfig <: BenchmarkConfig
    n_x::Int64
    n_u::Int64

    RandLinConfig(n_x::Int64, n_u::Int64)
        #Validate parameters
        @assert (n_x, n_u) .> 0 "Number of state variables, inputs and outputs must be larger than 0"

        return new(n_x, n_u)
    end
end

"""
External constructor for retrieving default RandLinConfig instances.
# Arguments
- `id`: string for identifying default instances, with possible values: `"A"`,`"B"`,`"C"`
- `n_x`: override for ``n_x``
- `n_u`: override for ``n_u``
"""
function RandLinConfig(id::String; n_x = nothing, n_u = nothing)
    params = matread(artifact"pH_RandLinConfig_"*id)    
    n_x == nothing ? n_x = params["n_x"] :
    n_u == nothing ? n_u = params["n_u"] :

    return RandLinConfig(n_x, n_u)
end

"""
Method for constructing the system matries.
# Arguments
- `config`: `RandLinConfig` instance
# Output
- `system`: Named tuple containing the system matrices in natural or pH form
"""
function construct_system(config::RandLinConfig; pH_form = false)
    function rand_SS(n) #Skew symmetric
        M = randn(n,n)
        return (M - M')/2
    end
    
    function rand_SPSD(n) #Symmetric positive semi-definite
        M = randn(n,n)
        return M * M'
    end

    function rand_SPD(n) #Symmetric positive definite
        M = rand_SPSD(n)
        return (M + M')/2
    end

    E = rand_SPD(config.n_x)
    Q = rand_SPD(config.n_x)

    Gamma = rand_SS(config.n_x + config.n_u)
    J = Gamma[1:n_x,1:n_x]
    G = Gamma[1:n_x,n_x+1:end]
    N = Gamma[n_x+1:end,n_x+1:end]

    W = rand_SPSD(config.n_x + config.n_u)
    R = W[1:n_x,1:n_x]
    P = W[1:n_x,n_x+1:end]
    S = W[n_x+1:end,n_x+1:end]

    return pH_form ? 
           (E = E, J = J, R = R, Q = Q, G = G, P = P, S = S, N = N) :
           (E = E, A = (J - R)*Q, B = G - P, C = (G + P)*Q, D = S + N)
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
- `Discretization`: detailed description of the discretization procedure and conversion to port-Hamiltonian form;
- `Interface`: section importing the docstrings from the corresponding `<Model>.jl` file;
- `References`: reference section in `BibTeX` format.

To contribute to the documentation, add a file based on the example below at `/docs/src/<Model>.md` and add `"<Model>.md"` to the `"Benchmark Systems"` list in `/docs/make.jl`. Please use the current 'Benchmark System' pages for stylistic reference.
````markdown
# Random Linear pH-System

## Description
This benchmark is a linear pH-system, of the following form:
```math
```
where the matrices are randomly generated as follows:
```math
``` 

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
