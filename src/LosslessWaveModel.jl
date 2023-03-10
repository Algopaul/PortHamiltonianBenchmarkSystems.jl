using LazyArtifacts

"""
    LosslessWaveModelConfig

Configuration for the lossless wave model benchmark system. Construct with `LosslessWaveModelConfigFactory`.

# Arguments
- `n`: System dimension
"""
struct LosslessWaveModelConfig <: BenchmarkConfig
    n::Int
    function LosslessWaveModelConfig(; n::Int)
        @assert n in [516741, 84225, 21701, 64013, 10841, 3189, 3825, 83521, 14089] "No data available for n"
        return new(n)
    end
end

"""
    LosslessWaveModelConfigFactory(; shape::String, size::String)

# Arguments
- `shape` (String): Shape of the domain. One of "disc", "L", "rectangle".
- `size` (String): Size of the domain. One of "small", "medium", "large".
"""
function LosslessWaveModelConfigFactory(; shape::String, size::String)
    n = losslesswavedimension(shape, size)
    return LosslessWaveModelConfig(; n)
end

function losslesswavedimension(shape::String, size::String)
    if shape == "disc" || shape == "Disc"
        return LosslessWaveDiscSizes[size]
    elseif shape == "L"
        return LosslessWaveLSizes[size]
    elseif shape == "rectangle" || shape == "Rectangle"
        return LosslessWaveRectSizes[size]
    else
        error("No data available for shape = $shape")
    end
end

LosslessWaveDiscSizes = Dict(["small" => 21701, "medium" => 84225, "large" => 516741])
LosslessWaveLSizes = Dict(["small" => 3189, "medium" => 10841, "large" => 64013])
LosslessWaveRectSizes = Dict(["small" => 3825, "medium" => 14089, "large" => 83521])

function construct_system(hmc::LosslessWaveModelConfig)
    data = load_losslesswave_data(hmc.n)
    return (
        M_p = data["M_p"],
        M_q = data["M_q"],
        M_b = data["M_b"],
        B = data["B"],
        D = data["D"],
    )
end

function load_losslesswave_data(n)
    if n == 21701
        data = artifact"Wave_Neumann_Disc_15266_6183_252"
        dd = matread(joinpath(data, "Wave_Neumann_Disc_15266_6183_252.mat"))
    elseif n == 516741
        data = artifact"Wave_Neumann_Disc_367930_147551_1260"
        dd = matread(joinpath(data, "Wave_Neumann_Disc_367930_147551_1260.mat"))
    elseif n == 84225
        data = artifact"Wave_Neumann_Disc_59692_24029_504"
        dd = matread(joinpath(data, "Wave_Neumann_Disc_59692_24029_504.mat"))
    elseif n == 3189
        data = artifact"Wave_Neumann_L_2162_903_124"
        dd = matread(joinpath(data, "Wave_Neumann_L_2162_903_124.mat"))
    elseif n == 64013
        data = artifact"Wave_Neumann_L_45162_18247_604"
        dd = matread(joinpath(data, "Wave_Neumann_L_45162_18247_604.mat"))
    elseif n == 10841
        data = artifact"Wave_Neumann_L_7520_3081_240"
        dd = matread(joinpath(data, "Wave_Neumann_L_7520_3081_240.mat"))
    elseif n == 3825
        data = artifact"Wave_Neumann_Rectangle_2620_1085_120"
        dd = matread(joinpath(data, "Wave_Neumann_Rectangle_2620_1085_120.mat"))
    elseif n == 83521
        data = artifact"Wave_Neumann_Rectangle_59100_23821_600"
        dd = matread(joinpath(data, "Wave_Neumann_Rectangle_59100_23821_600.mat"))
    elseif n == 14089
        data = artifact"Wave_Neumann_Rectangle_9840_4009_240"
        dd = matread(joinpath(data, "Wave_Neumann_Rectangle_9840_4009_240.mat"))
    else
        error("No data available for n = $n")
    end
    return dd
end

export LosslessWaveModelConfig, LosslessWaveModelConfigFactory
