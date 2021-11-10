module PortHamiltonianBenchmarkSystems

using LinearAlgebra, SparseArrays, HDF5
"""
gugercin_pH_msd_chain(; n_cells=50::Int, m=2::Int, c_i=1.0, m_i=4.0, k_i=4.0)
This function returns the port Hamiltonian mass-spring-damper system described in
S. Gugercin et al.:
      Structure-preserving tangential interpolation for model reduction of
      port-Hamiltonian systems
# Arguments
- n_cells describes the number of masses. The system dimension is 2n_cells
- c_i stears the amount of damping
- m_i is the weight of the masses
- k_i is the stiffness of the springs
"""
function gugercin_pH_msd_chain(;
    n_cells=50::Int,
    m=2::Int,
    c_i=1.0,
    m_i=4.0,
    k_i=4.0
  )
  n=2*n_cells;
  # B is initialized as dense matrix. Since all results of transfer function
  # computations will lead to dense results.
  B=zeros(n,m);
  [B[2*i,i]=1.0 for i in 1:m]
  J=spzeros(n,n);
  [J[i,i+1]=1.0 for i in 1:2:(n-1)];
  J=J-J'
  # Set constants.
  R=spzeros(n,n);
  [R[i,i]=c_i for i in 2:2:n]
  Q=spzeros(n,n);
  Q[1,1]=k_i
  [Q[i,i]=2*k_i for i in 3:2:(n-1)]
  [Q[i,i]=1/m_i for i in 2:2:n]
  [Q[i,i+2]=-k_i for i in 1:2:(n-2)]
  [Q[i+2,i]=-k_i for i in 1:2:(n-2)]
  return J, R, Q, B
end

end # module
