# Contribution

## Modus Operandi
This benchmark collection is driven by the active support of the port-Hamiltonian community. If your research has lead to port-Hamiltonian models that may be relevant for this collection, we would be happy to include them. If you wish to contribute to the project directly, please consult this page, fork our [Git repository](https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl/) and issue a pull request. Alternatively, feel free to contact us via [e-mail](mailto:schwerdt@math.tu-berlin.de) to discuss our potential collaboration. We are happy to receive reference implementations in other languages and reimplement them in Julia for this package.

All entries in the collection consist of code for generating the system matrices, a comprehensive set of tests and a documentation page. Exact instructions for each of these parts are provided in the sections below. Some of the design choices may still be subject to change, but efforts are currently being made to settle on a final structure. Any subsequently required changes to community contributions will be handled by us.

## Code
This package consists of a single module, containing the following five elements:
1. `PHSystem`: parametric composite type for storing system matrices in standard port-Hamiltonian form,
2. `<System>Config`: (parametric) composite type for storing parameter sets for some system,
3. `<System>Config(id::String)`: external constructor returing default `<System>Config` instances, based on some `id`,
4. `construct_system(config::<System>Config)`: method returning system matrices in "natural" form, based on some `<System>Config` instance,
5. `PHSystem(config::<System>Config)`: external constructor returning `PHSystem` instances, based on some `<System>Config` instance.
The last four elements are repeated for each benchmark system and are stored together in `/src/<System>.jl`, along with their respective docstrings in Markdown format. All `<System>.jl` files are included in `/src/PortHamiltonianBenchmarkSystems.jl`.

To contribute to the code, simply add a file based on the example below at `/src/<System>.jl` and add a corresponding `include` statement to `/src/PortHamiltonianBenchmarkSystems.jl`. Some best practices for more complex systems are also given below.
```@raw html
<details><summary>Code Example</summary>
```
```julia
export RandLinConfig

"""
Composite type describing a linear port-Hamiltonian system, where matrices ``J,\ R,\ G`` are random dense matrices of the correct structure (R positive semi-definite, J, G skew symmetric) with mean 0 and variance 1.
# Arguments
- n_x: number of state variables
- n_p: number of input and output ports
"""
struct RandLinConfig
    n_x::UInt64
    n_p::UInt64
end

"""
Method for constructing the system matries in "natural" form.
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
    G = (M3 - M3')/sqrt(2)

    return (J = J, R = R, G = G)
end

function PHSystem(config::RandLinConfig)
    J, R, G = construct_system(config)

    E = sparse(1.0I,size(J)...)
    Q = sparse(1.0I,size(J)...)
    P = spzeros(size(G)...)
    N = spzeros(config.n_p, config.n_p)
    S = spzeros(size(N)...)

    return PHSystem(E, J, R, Q, G, P, S, N)
end
```
```@raw html
</details>
```

```@raw html
<details><summary>Best Practices</summary>
```
While it may seem unnecessary in a simple case like this, it is important to include an internal constructor in `<System>Config`, to perform input validation:
```julia
struct RandLinConfig
    n_x::UInt64
    n_p::UInt64

    RandLinConfig(n_x::UInt64, n_p::UInt64)
        @assert n_x > 0 "Number of state variables must be larger than 0"
        @assert n_p > 0 "Number of ports must be larger than 0"
        return new(n_x, n_p)
    end
end
```
Similarly, while it does not make much sense to save default configurations for this random system, in reality one should provide a few default parameter sets. Once some adequate parameter sets have been determined and implemented in the form shown below, we will store them as MAT-files on our [Zenodo](https://zenodo.org/) and update the code to retrieve them as Julia [artifacts](https://docs.julialang.org/en/v1/stdlib/Artifacts/) when needed.
```julia
"""
Constructor for retrieving default `RandLinConfig` instances.
# Arguments
- `id`: identifier for default parameter set
# Output
- `system`: `RandLinConfig` instance
"""
function RandLinConfig(id::String)
    Random.seed!(0)

    if id == "A"
        return RandLinConfig(5, 3)
    elseif id == "B"
        return RandLinConfig(6, 2)
    elseif id == "C"
        return RandLinConfig(7, 4)
    else
        error("Invalid id!")
    end
end
```
```@raw html
</details>
```

## Tests
To guarantee that all merged code is in a working state, automatic test pipelines have been set up for this repository. The test are run through GitHub Actions, using the Julia [Test](https://docs.julialang.org/en/v1/stdlib/Test/) module. 

Every benchmark model has its own `/test/<System>Tests.jl` file, which is included in `/test/runtests.jl`. While the goal is to achieve near total code coverage, it is not necessary to 'test everything'. We suggest to at least test the following for a variety of configurations:
- System matrix sizes,
- System transfer functions,
- Correspondence between 'natural' and port-Hamiltonian system matrices.

To contribute to the tests, add a file based on the example below at `/test/<System>Tests.jl` and add a corresponding `include` statement to `/test/runtests.jl`.
```@raw html
<details><summary>Test Example</summary>
```
```julia
@testset "RandLin" begin
    ids = ["A", "B", "C"]

    for id in ids
        config = RandLinConfig(id)
        J, R, G = construct_system(config)
        sys = PHSystem(config)

        #System matrix sizes
        @test size(J) == (n_x, n_x)
        @test size(R) == (n_x, n_x)
        @test size(G) == (n_x, n_p)

        #Correspondence between "natural" and pH form
        @test sys.E == sparse(1.0I,size(J)...)
        @test (sys.S, sys.N) .== spzeros(config.n_p, config.n_p)
        @test (sys.J - sys.R) * sys.Q == J - R
        @test (sys.G - sys.P) == G
        @test (sys.G + sys.P)' == G'
    end
end
```
```@raw html
</details>
```

## Documentation
The documentation for this package is built using [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl). The `/docs/make.jl` script uses the Markdown files in `/docs/src` and the images in `/docs/src/assets` to build a documentation webpage in `/docs/build`. Equations are rendered using [``\KaTeX``](https://katex.org/), which is invoked in Markdown by the `math` environment. The webpage can be loaded locally by running the commands below in the Julia REPL. The documentation should then be accessible from the returned `localhost` port.
```julia
using LiveServer
include("/docs/make.jl")
serve(dir="/docs/build")
```
Each benchmark model is documented in a separate file, containing the following sections:
- Description: mathematical description of the model,
- Derivation: detailed description of the conversion from mathematical model to port-Hamiltonian system,
- Interface: section importing the docstrings from the corresponding `<System>.jl` file,
- References: reference section in [BibTeX](http://www.bibtex.org/) format.
To contribute to the documentation, add a file based on the example below at `/docs/src/<System>.md` and add `"<System>.md"` to the `"Benchmark Systems"` list in `/docs/make.jl`. Please refer to the current "Benchmark Systems" pages for stylistic reference.
```@raw html
<details><summary>Documentation Example</summary>
```
````markdown
# Random Linear pH-System

## Description
This benchmark is a linear port-Hamiltonian system of the form:
```math
\begin{align*}
    E\dot{x} &= (J-R)Qx + (G-P)u,\\
    y &= (G+P)^HQx + (S+N)u,
\end{align*}
```
where ``E,Q=I``, ``P,S,N=0`` and ``J,\ R,\ G`` are random dense matrices of the correct structure (R positive semi-definite, J, G skew symmetric) with mean 0 and variance 1 [Sabbadini2022](#References).

## Derivation
This system discrete, random, and in port-Hamiltonian form automatically.

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
@article{Sabbadini2022,
  title = "On random port-Hamiltonian Systems",
  journal = "We publish anything"
}
```
````
```@raw html
</details>
```