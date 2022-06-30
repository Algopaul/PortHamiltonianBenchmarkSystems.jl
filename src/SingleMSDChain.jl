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
struct SingleMSDConfig{TC, TM, TK} <: BenchmarkConfig
  n_cells::Int
  io_dim::Int
  c::TC
  m::TM
  k::TK
  function SingleMSDConfig(n_cells::Int, io_dim::Int, c::TC, m::TM, k::TK) where {TC, TM, TK}
    @assert n_cells > 0 "number of cells must be positive"
    @assert io_dim > 0 "number of inputs and outputs must be positive"
    @assert c >= 0 "damping cannot be negative"
    @assert m >= 0 "masses cannot be negative"
    @assert k >= 0 "stiffnesses cannot be negative"
    return new{TC, TM, TK}(n_cells, io_dim, c, m, k)
  end
end

function SingleMSDConfig()
  return SingleMSDConfig(100, 2, 1.0, 4.0, 4.0)
end

function construct_system(config::SingleMSDConfig{TC, TM, TK}) where {TC <: Number, TM <: Number, TK <: Number}
  @unpack n_cells, io_dim, c, m, k = config
  n=2*n_cells;
  # B is initialized as dense matrix. Since all results of transfer function
  # computations will lead to dense results.
  B=zeros(n, io_dim);
  [B[2*i,i]=1.0 for i in 1:io_dim]
  J=spzeros(n,n);
  [J[i,i+1]=1.0 for i in 1:2:(n-1)];
  J=J-J'
  # Set constants.
  R=spzeros(n,n);
  [R[i,i]=c for i in 2:2:n]
  Q=spzeros(n,n);
  Q[1,1]=k
  [Q[i,i]=2*k for i in 3:2:(n-1)]
  [Q[i,i]=1/m for i in 2:2:n]
  [Q[i,i+2]=-k for i in 1:2:(n-2)]
  [Q[i+2,i]=-k for i in 1:2:(n-2)]
  return (J=J, R=R, Q=Q, B=B)
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
function gugercin_pH_msd_chain(;
    n_cells=50,
    m=2,
    c_i=1.0,
    m_i=4.0,
    k_i = 4.0
  )
  config = SingleMSDConfig(n_cells, m, c_i, m_i, k_i)
  return construct_system(config)
end

export SingleMSDConfig, gugercin_pH_msd_chain
