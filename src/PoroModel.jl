using MAT
using LazyArtifacts

"""
This struct configures port Hamiltonian poroelasticity systems described in
    Altmann, Mehrmann, Unger: Port-Hamiltonian Formulations of Poroelastic
    Network Models
# Arguments
- `n`: System dimension (can only be either: 320, 980, or 1805). Default = 980.
- `m`: Number of inputs/outputs (can only be 1 or 2). Default = 1.
- `rho`: density. Default = `1e-3`.
- `alpha`: Biot-Willis fluid-solid coupling coefficient. Default = 0.79.
- `bm`: Biot-Modulus. Default = `1/7.8e3`.
- `kappanu`: Quotient `kappa/Nu`, where `kappa` denotes the permeability and `nu` denotes the fluid viscosity. Default = 633.33.
- `eta`: artificial damping coefficient. Default = `1e-4`.
"""
struct PoroElasticityConfig <: BenchmarkConfig
    n::Int
    m::Int
    rho::Float64
    alpha::Float64
    bm::Float64
    kappanu::Float64
    eta::Float64
    function PoroElasticityConfig(;
        n::Int = 980,
        m::Int = 1,
        rho::Float64 = 1e-3,
        alpha::Float64 = 0.79,
        bm::Float64 = 1 / 7.80e3,
        kappanu::Float64 = 633.33,
        eta::Float64 = 1e-4,
    ) where {}
        @assert n in [320, 980, 1805] "n must be one of 320, 980, or 1805"
        @assert m in [1, 2] "m must be 1 or 2 for poroelasticity model"
        @assert rho > 0 "rho must be positive"
        @assert alpha > 0 "alpha must be positive"
        @assert bm > 0 "bm must be positive"
        @assert kappanu > 0 "kappanu must be positive"
        return new(n, m, rho, alpha, bm, kappanu, eta)
    end
end

function construct_system(config::PoroElasticityConfig)
    (; n, m, rho, alpha, bm, kappanu, eta) = config
    Y, D, M, K, Bp, Bf, A = load_poro_raw_data(n = n)
    Y = rho * sparse(Y)
    D = alpha * sparse(D)
    M = 1 / bm .* sparse(M)
    K = kappanu * sparse(K)
    A = sparse(A)
    Bp = Bp'
    Bf = Bf'
    n = size(A, 1)
    
    l = size(M, 1)
    E = [Y spzeros(n, n + l); spzeros(n, n) A spzeros(n, l); spzeros(l, n + n) M]
    J = [spzeros(n, n) -A D'; A spzeros(n, n + l); -D spzeros(l, n + l)]
    R = [spzeros(n, 2 * n + l); spzeros(n, 2 * n + l); spzeros(l, 2 * n) K] + eta * I
    
    if m == 1
        G = [zeros(n, 1); Bf; Bp]
    elseif m == 2
        G = [Bf zeros(n,1); zeros(n,2); zeros(l,1) Bp];
    else
        error("m must be 1 or 2 for poroelasticity model")
    end

    return (E = E, J = J, R = R, G = G)
end

function load_poro_raw_data(; n = 980)
    poro_data = artifact"ph-poromodelsbasedata"
    matfile = joinpath(poro_data, "PH-PoroModelsBaseData/poro-n$n.mat")
    dd = matread(matfile)
    return dd["Y"], dd["D"], dd["M"], dd["K"], dd["Bp"], dd["Bf"], dd["A"]
end

function PHSystem(config::PoroElasticityConfig)
    E, J, R, G = construct_system(config)
    n, m = size(G)
    Q = I(n)
    P = spzeros(n, m)
    S = spzeros(m, m)
    N = spzeros(m, m)
    return PHSystem(E, J, R, Q, G, P, S, N)
end

"""
    poro_elasticity_model(; n = 980, m = 1, rho = 1e-3, alpha = 0.79, bm = 1/7.80e3, kappanu = 633.33, eta = 1e-4)

This function returns a port-Hamiltonian model of linear poroelasticity in a
bounded Lipschitz domain as described in
    Altmann, Mehrmann, Unger: Port-Hamiltonian Formulations of Poroelastic
    Network Models
# Arguments
- `n`: System dimension (can only be either: 320, 980, or 1805). Default = 980.
- `m`: Number of inputs/outputs (can only be 1 or 2). Default = 1.
- `rho`: density. Default = `1e-3`.
- `alpha`: Biot-Willis fluid-solid coupling coefficient. Default = 0.79.
- `bm`: Biot-Modulus. Default = `1/7.8e3`.
- `kappanu`: Quotient `kappa/Nu`, where `kappa` denotes the permeability and `nu` denotes the fluid viscosity. Default = 633.33.
- `eta`: artificial damping coefficient. Default = `1e-4`.
# Outputs
- ``E, J, R, B``, matrices to construct the transfer function ``H(s) = B^\\mathsf{T}(sE-(J-R))^{-1}B)``
"""
function poro_elasticity_model(;
    n = 980,
    m = 1,
    rho = 1e-3,
    alpha = 0.79,
    bm = 1 / 7.80e3,
    kappanu = 633.33,
    eta = 1e-4,
)
    config = PoroElasticityConfig(n=n, m=m, rho=rho, alpha=alpha, bm=bm, kappanu=kappanu, eta=eta)
    return construct_system(config)
end

export PoroElasticityConfig, poro_elasticity_model
