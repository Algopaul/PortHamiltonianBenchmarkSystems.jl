using UnPack

"""
This struct configures port Hamiltonian mass-spring-damper systems described in
S. Gugercin et al.:
      Structure-preserving tangential interpolation for model reduction of
      port-Hamiltonian systems
# Arguments
- `n_cells`: The number of masses. The system dimension is `2n_cells`
- `c_i`: The amount of damping
- `m_i`: The weight of the masses
- `k_i`: The stiffness of the springs
# Outputs
Matrices: ``J, R, Q, B``. The resulting transfer function is ``H(s) = B^\\mathsf{T} Q  (sI-(J-R)Q)^{-1}B``.
"""
struct SingleMSDConfig{TC,TM,TK} <: BenchmarkConfig
    n_cells::Int
    io_dim::Int
    c::TC
    m::TM
    k::TK
    function SingleMSDConfig(
        n_cells::Int=50,
        io_dim::Int=2,
        c::TC=1.0,
        m::TM=4.0,
        k::TK=4.0,
    ) where {TC,TM,TK}
        @assert n_cells > 0 "number of cells must be positive"
        @assert io_dim > 0 "number of inputs and outputs must be positive"
        @assert io_dim <= n_cells "number of inputs and outputs must be less than or equal to the number of cells"
        msd_check_constants(c, "damping")
        msd_check_constants(m, "masses")
        msd_check_constants(k, "stiffness")
        return new{TC,TM,TK}(n_cells, io_dim, c, m, k)
    end
end

function msd_check_constants(c::Number, name)
    @assert c >= 0 "$name cannot be negative"
end

function msd_check_constants(c::AbstractVector, name)
    @assert all(c .>= 0) "$name cannot have any negative entries"
end

function construct_system(config::SingleMSDConfig)
    @unpack n_cells, io_dim, c, m, k = config
    n = 2 * n_cells
    # B is initialized as dense matrix. Since all results of transfer function
    # computations will lead to dense results.
    B = zeros(n, io_dim)
    [B[2 * i, i] = 1.0 for i = 1:io_dim]
    J = spzeros(n, n)
    [J[i, i + 1] = 1.0 for i = 1:2:(n - 1)]
    J = J - J'
    # Set constants.
    R = msd_construct_R(n_cells, c)
    Q = msd_construct_Q(n_cells, k, m)
    return (J = J, R = R, Q = Q, B = B)
end

function msd_construct_R(n_cells, c)
    n = 2n_cells
    R = spzeros(n, n)
    for (j, i) in enumerate(2:2:n)
        R[i, i] = msd_vecint(c, j)
    end
    return R
end

function msd_construct_Q(n_cells, k, m)
    n = 2n_cells
    Q = spzeros(n, n)
    n = size(Q, 1)
    for (j, i) in enumerate(1:2:(n - 3))
        msd_construct_Q_add_k_stencil(Q, msd_vecint(k, j), i)
    end
    Q[end - 1, end - 1] += msd_vecint(k, n_cells)
    for (j, i) in enumerate(2:2:n)
        Q[i, i] = 1 / msd_vecint(m, j)
    end
    return Q
end

function msd_vecint(x::Number, ::Any)
    return x
end

function msd_vecint(x::AbstractVector, i)
    return x[i]
end

function msd_construct_Q_add_k_stencil(Q, k, i)
    @views Q[i:(i + 2), i:(i + 2)] .+= [k 0 -k; 0 0 0; -k 0 k]
end

function PHSystem(config::SingleMSDConfig)
    J, R, Q, G = construct_system(config)
    n, m = size(G)
    E = I(n)
    P = spzeros(n, m)
    S = spzeros(m, m)
    N = spzeros(m, m)
    return PHSystem(E, J, R, Q, G, P, S, N)
end

"""
    gugercin_pH_msd_chain(; n_cells=50, m=2, c_i=1.0, m_i=4.0, k_i=4.0)

This function returns the port Hamiltonian mass-spring-damper system described in
S. Gugercin et al.:
      Structure-preserving tangential interpolation for model reduction of
      port-Hamiltonian systems
# Arguments
- `n_cells`: The number of masses. The system dimension is `2n_cells`
- `c_i`: The amount of damping
- `m_i`: The weight of the masses
- `k_i`: The stiffness of the springs
# Outputs
Matrices: ``J, R, Q, B``. The resulting transfer function is ``H(s) = B^\\mathsf{T} Q  (sI-(J-R)Q)^{-1}B``.
"""
function gugercin_pH_msd_chain(; n_cells = 50, m = 2, c_i = 1.0, m_i = 4.0, k_i = 4.0)
    @warn """Deprecated, use SingleMSDConfig("Gugercin") instead"""
    config = SingleMSDConfig(n_cells, m, c_i, m_i, k_i)
    return construct_system(config)
end

"""
    SingleMSDConfig(id::String)

External constructor providing various default instances of SingleMSDConfig.
# Arguments
- `id`: The identifier of the desired configuration.
# Outputs
- `config`: Instance of `SingleMSDConfig`.
"""
function SingleMSDConfig(id::String)
    if id == "Gugercin"
        return SingleMSDConfig(50, 2, 1.0, 4.0, 4.0)
    else
        error("Unknown benchmark id: " + id)
    end
end

export SingleMSDConfig, gugercin_pH_msd_chain
