using MAT
using SparseArrays
using LazyArtifacts

export Elasticity2DAFWConfig

"""
Composite type describing port-Hamiltonian elasticity systems described in A. 
Brugnoli: A port-Hamiltonian formulation of flexible structures. Modelling 
and structure-preserving finite element discretization.
# Arguments
- `n`: System dimension, one of ``\\{1260, 1880, 4920\\}``.
"""
struct Elasticity2DAFWConfig <: BenchmarkConfig
    n::Int
    function Elasticity2DAFWConfig(n::Int)
        @assert n in [1260, 1880, 4920] "n must be one of 1260, 1880, or 4920."
        return new(n)
    end
end

"""
External constructor providing the default instance of Elasticity2DAFWConfig,
namely n = 1880.
"""
function Elasticity2DAFWConfig()
    return Elasticity2DAFWConfig(1880)
end

function construct_system(config::Elasticity2DAFWConfig)
    E, J, B, coord_u = load_el2Dafw_raw_data(config.n)
    return (E = E, J = J, G = B)
end

function load_el2Dafw_raw_data(n)
    poro_data = artifact"elasticity_model"
    matfile = joinpath(poro_data, "el2Dafw-n$n.mat")
    dd = matread(matfile)
    return dd["E"], dd["J"], dd["B"], dd["x_u"]
end

function PHSystem(config::Elasticity2DAFWConfig)
    E, J, G = construct_system(config)
    n, m = size(G)
    R = spzeros(n, n)
    Q = sparse(1.0I, size(J)...)
    P = spzeros(n, m)
    S = spzeros(m, m)
    N = spzeros(m, m)
    return PHSystem(E, J, R, Q, G, P, S, N)
end