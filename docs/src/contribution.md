# Contribution
## Code
This package consists of a single module, containing the following six elements:
- `PHSystem`: parametric composite type for storing system matrices in port-Hamiltonian form, with a single input validating internal constructor;
- `BenchmarkConfig`: abstract type for grouping all `<Model>Config` types;
- `<Model>Config <: BenchmarkConfig`: (parametric) composite type containing all parameters for some model, with a single input validating internal constructor;
- `<Model>Config`: external constructor providing several default `<Model>Config` instances based on some identifier
- `construct_system`: method for constructing system matrices in 'natural' form based on some `<Model>Config` instance;
- `PHSystem`: external constructor for constructing system matrices in port-Hamltonian form.
The last four elements are repeated for each benchmark model and are stored together in `src/<Model>.jl`. All `<Model>.jl` files are included in `src/PortHamiltonianBenchmarkSystems.jl`, which also contains the type declarations for `PHSystem` and `BenchmarkConfig`.

To contribute to the code, simply add a file based on the template below at `/src/<Model>.jl` and add a corresponding `include` statement to `src/PortHamiltonianBenchmarkSystems.jl`. The docstrings should be written in `Markdown` format, as shown below.
```julia
export <Model>Config

"""
<description>
# Arguments
- <args>: <types>, <descriptions>
"""
struct <Model>Config <: BenchmarkConfig
    <parameters>::<types>

    <Model>Config(<parameters>::<types>)
        #Validate <parameters>

        return new(<parameters>)
    end
end

"""
External constructor providing various default <Model> configurations.
# Arguments
- `id`: <type> to identify a default configurations, with possible values: <values>
"""
function <Model>Config(id::<type>)
    #Fetch <parameters> based on id

    return <Model>Config(<parameters>)
end


"""
Method for constructing the 'natural' system matries.
# Arguments
- `config`: `<Model>Config` instance
# Output
- `system`: Named tuple containing sparse matrices `<names>`
"""
function construct_system(config::<Model>Config)
    #Construct <natural system matrices> based on config

    return (<names> = <natural system matrices>) #Named tuple
end

function PHSystem(config::<Model>Config)
    <natural system matrices> = construct_system(config)

    #Construct <pH system matrices> based on <natural system matrices>
    
    return PHSystem(<pH system matrices>)
end
```
## Tests
To guarantee that all merged code is in a working state, automatic test pipelines have been set up for this repository. The test are run through GitHub Actions, using the `Test.jl` package. 

Every benchmark model has its own `/test/<Model>Tests.jl` file, which is included in `/test/runtests.jl`. While the goal is to achieve near total code coverage, it is not necessary to 'test everything'. We suggest to at least test the following for a variety of configurations:
- System matrix sizes;
- System transfer functions;
- Correspondence between 'natural' and port-Hamiltonian system matrices.

To contribute to the tests, add a file based on the template below at `/test/<Model>Tests.jl` and add a corresponding `include` statement to `/test/runtests.jl`. Please refer to the current tests for more detailed examples.
```julia
@testset "<Model>" begin
    @test #Test 1
    @test #Test ...
    @test #Test N
end
```
## Documentation
The documentation for this package is built using `Documenter.jl`. The `/docs/make.jl` script uses the `Markdown` files in `/docs/src` and the images in `/docs/src/assets`to build a documentation webpage in `/docs/build`. The webpage can be loaded locally by running `make.jl` and then `LiveServer.serve(dir=/docs/build`). The documentation should then be accessible from the returned `http://localhost` port.

Each benchmark model is documented in a separate file, containing the following sections:
- `Description`: mathematical description of the model;
- `Discretization`: detailed description of the discretization procedure and conversion to port-Hamiltonian form;
- `Interface`: section importing the docstrings from the corresponding `<Model>.jl` file;
- `References`: reference section in `BibTeX` format.

To contribute to the documentation, add a file based on the template below at `/docs/src/<Model>.md` and add `"<Model>.md"` to the `"Benchmark Systems"` list in `/docs/make.jl`. Please use the current 'Benchmark System' pages for stylistic reference.
````markdown
# <Model>

## Description

## Discretization

## Interface
```@docs
<Model>Config
```
```@docs
<Model>Config(id::<type>)
```
```@docs
construct_system(config::<Model>Config)
```

## References
```LaTeX

```
````