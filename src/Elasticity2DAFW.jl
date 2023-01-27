using MAT
using SparseArrays
using LazyArtifacts

"""
This struct configures port Hamiltonian elasticity systems described in
    A. Brugnoli: A port-Hamiltonian formulation of flexible structures.
    Modelling and structure-preserving finite element discretization
# Arguments
- `n`: System dimension (can only be either: 1260, 1880, or 4920). Default = 1880.
"""
struct Elasticity2DAFWConfig <: BenchmarkConfig
    n::Int
    function Elasticity2DAFWConfig(; n::Int = 1260) where {}
        @assert n in [1260, 1880, 4920] "n must be one of 1260, 1880, or 4920."
        return new(n)
    end
end

function construct_system(config::Elasticity2DAFWConfig)
    (; n) = config
    E, J, B, coord_u = load_el2Dafw_raw_data(n = n)

    n = size(E, 1)
    m = size(B, 2)
    R = spzeros(n, n)
    return (E = E, J = J, R = R, G = B)
end

function load_el2Dafw_raw_data(; n = 1260)
    poro_data = artifact"elasticity_model"
    matfile = joinpath(poro_data, "el2Dafw-n$n.mat")
    dd = matread(matfile)
    return dd["E"], dd["J"], dd["B"], dd["x_u"]
end

function PHSystem(config::Elasticity2DAFWConfig)
    E, J, R, G = construct_system(config)
    n, m = size(G)
    Q = I(n)
    P = spzeros(n, m)
    S = spzeros(m, m)
    N = spzeros(m, m)
    return PHSystem(E, J, R, Q, G, P, S, N)
end

export Elasticity2DAFWConfig
