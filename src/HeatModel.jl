using LazyArtifacts

"""
    HeatModelConfig

Configuration for the heat model benchmark system. Construct with `HeatModelConfigFactory`.

# Arguments
- `n`: System dimension
"""
struct HeatModelConfig <: BenchmarkConfig
    n::Int
    function HeatModelConfig(; n::Int)
        @assert n in [516741, 84225, 21701, 64013, 10841, 3189, 3825, 83521, 14089] "No data available for n"
        return new(n)
    end
end

"""
    HeatModelConfigFactory(; shape::String, size::String)

# Arguments
- `shape` (String): Shape of the domain. One of "disc", "L", "rectangle".
- `size` (String): Size of the domain. One of "small", "medium", "large".
"""
function HeatModelConfigFactory(; shape::String, size::String)
    n = heatdimension(shape, size)
    return HeatModelConfig(; n)
end

function heatdimension(shape::String, size::String)
    if shape == "disc" || shape == "Disc"
        return HeatDiscSizes[size]
    elseif shape == "L"
        return HeatLSizes[size]
    elseif shape == "rectangle" || shape == "Rectangle"
        return HeatRectSizes[size]
    else
        error("No data available for shape = $shape")
    end
end

HeatDiscSizes = Dict(["small" => 21701, "medium" => 84225, "large" => 516741])
HeatLSizes = Dict(["small" => 3189, "medium" => 10841, "large" => 64013])
HeatRectSizes = Dict(["small" => 3825, "medium" => 14089, "large" => 83521])

function construct_system(hmc::HeatModelConfig)
    data = load_heat_data(hmc.n)
    return (
        M_T = data["M_T"],
        M_Q = data["M_Q"],
        M_b = data["M_b"],
        D = data["D"],
        J = data["J"],
    )
end

function load_heat_data(n)
    if n == 516741
        data = artifact"Heat_Neumann_Disc_147551_367930_1260"
        dd = matread(joinpath(data, "Heat_Neumann_Disc_147551_367930_1260.mat"))
    elseif n == 84225
        data = artifact"Heat_Neumann_Disc_24029_59692_504"
        dd = matread(joinpath(data, "Heat_Neumann_Disc_24029_59692_504.mat"))
    elseif n == 21701
        data = artifact"Heat_Neumann_Disc_6183_15266_252"
        dd = matread(joinpath(data, "Heat_Neumann_Disc_6183_15266_252.mat"))
    elseif n == 64013
        data = artifact"Heat_Neumann_L_18247_45162_604"
        dd = matread(joinpath(data, "Heat_Neumann_L_18247_45162_604.mat"))
    elseif n == 10841
        data = artifact"Heat_Neumann_L_3081_7520_240"
        dd = matread(joinpath(data, "Heat_Neumann_L_3081_7520_240.mat"))
    elseif n == 3189
        data = artifact"Heat_Neumann_L_903_2162_124"
        dd = matread(joinpath(data, "Heat_Neumann_L_903_2162_124.mat"))
    elseif n == 3825
        data = artifact"Heat_Neumann_Rectangle_1085_2620_120"
        dd = matread(joinpath(data, "Heat_Neumann_Rectangle_1085_2620_120.mat"))
    elseif n == 83521
        data = artifact"Heat_Neumann_Rectangle_23821_59100_600"
        dd = matread(joinpath(data, "Heat_Neumann_Rectangle_23821_59100_600.mat"))
    elseif n == 14089
        data = artifact"Heat_Neumann_Rectangle_4009_9840_240"
        dd = matread(joinpath(data, "Heat_Neumann_Rectangle_4009_9840_240.mat"))
    else
        error("No data available for n = $n")
    end
    return dd
end

export HeatModelConfig, HeatModelConfigFactory
