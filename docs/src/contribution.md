# Contribution
## Code
This package consists of a module, containing the following six elements:
- `PHSystem`: parametric composite type for storing system matrices in port-Hamiltonian form, with a single input validating internal constructor
- `BenchmarkConfig`: abstract type for grouping all `ModelConfig` types
- `ModelConfig <: BenchmarkConfig`: (parametric) composite type containing all parameters for some model, with a single input validating internal constructor
- `ModelConfig`: External constructor providing several default `ModelConfig` instances based on some identifier
- `construct_system`: Method for constructing system matrices in 'natural' form based on some `ModelConfig` instance
- `PHSystem`: External constructor for constructing system matrices in port-Hamltonian form
The last four elements are repeated for each benchmark model and are stored together in `Model.jl`. All `Model.jl` files are included in `PortHamiltonianBenchmarkSystems.jl`, which also contains the type declarations for `PHSystem` and `BenchmarkConfig` along with some common package imports.

To contribute, simply add a `Model.jl` file based on the template shown below to `/src` and add a corresponding `include` statement to `src/PortHamiltonianBenchmarkSystems.jl`.
```julia
export ModelConfig

struct ModelConfig <: BenchmarkConfig
    <parameters>::<types>

    ModelConfig(<parameters>::<types>)
        #Validate <parameters>

        return new(<parameters>)
    end
end

function ModelConfig(id::<type>)
    #Fetch <parameters> based on id

    return ModelConfig(<parameters>)
end

function construct_system(config::ModelConfig)
    #Construct <natural system matrices> based on config

    return (<names> = <natural system matrices>) #Named tuple
end

function PHSystem(config::ModelConfig)
    <natural system matrices> = construct_system(config)

    #Construct <pH system matrices> based on <natural system matrices>
    
    return PHSystem(<pH system matrices>)
end
```
## Documentation