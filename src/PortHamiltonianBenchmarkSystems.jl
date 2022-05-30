module PortHamiltonianBenchmarkSystems

using LinearAlgebra, SparseArrays

include("IOFormats.jl")
include("Downloads.jl")

"""
   `gugercin_pH_msd_chain(; n_cells=50, m=2, c_i=1.0, m_i=4.0, k_i=4.0)`

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

"""
`poro_elasticity_model(;
    n = 980,
    rho = 1e-3,
    alpha = 0.79,
    M = 1/7.80e3,
    kappanu = 633.33,
    eta = 1e-4,
    force_download = false
  )`

This function returns a port-Hamiltonian model of linear poroelasticity in a
bounded Lipschitz domain as described in
    Altmann, Mehrmann, Unger: Port-Hamiltonian Formulations of Poroelastic
    Network Models
# Arguments
- `n`: System dimension (can only be either: 320, 980, or 1805). Default = 980.
- `rho`: density. Default = `1e-3`.
- `alpha`: Biot-Willis fluid-solid coupling coefficient. Default = 0.79.
- `bm`: Biot-Modulus. Default = `1/7.8e3`.
- `kappanu`: Quotient kappa/Nu, where kappa denotes the permeability and nu denotes the fluid viscosity. Default = 633.33.
- `eta`: artificial damping coefficient. Default = `1e-4`.
"""
function poro_elasticity_model(;
    n = 980,
    rho = 1e-3,
    alpha = 0.79,
    bm = 1/7.80e3,
    kappanu = 633.33,
    eta = 1e-4,
    force_download = false
  )
  Y, D, M, K, Bp, Bf, A = load_poro_raw_data(
    n=n,
    force_download = force_download
  )
  Y = rho*sparse(Y)
  D = alpha*sparse(D)
  M = 1/bm .* sparse(M)
  K = kappanu*sparse(K)
  A = sparse(A)
  Bp = Bp'
  Bf = Bf'
  n = size(A, 1);
  m = size(M, 1);
  E = [Y spzeros(n,n+m); spzeros(n,n) A spzeros(n,m); spzeros(m,n+n) M];
  J = [spzeros(n,n) -A D';A spzeros(n,n+m); -D spzeros(m,n+m)];
  R = [spzeros(n,2*n+m); spzeros(n,2*n+m); spzeros(m,2*n) K] + eta*I
  B = [zeros(n,1); Bf; Bp];
  return E, J, R, B
end

function load_poro_raw_data(;
    n = 980,
    force_download = false
  )
  filename = "poro-n$n.mat"
  url = "https://zenodo.org/record/5702554/files/poro-n$n.mat?download=1"
  if n == 980
    md5_hash = hex(0x2961a189be7049ffe2d476b18cb1f678)
  elseif n == 320
    md5_hash = hex(0x97afe8c34f0e9a56bbe86d0a51b7b626)
  elseif n == 1805
    md5_hash = hex(0xc61f6687da9cd26cbf2d880d7d3a9ac9)
  else
    throw(ArgumentError("Model size is either 320, 980, or 1805"))
  end
  download_system_data_if_required(
    filename,
    url,
    md5_hash,
    force_download = force_download
  )
  dd = loadMAT(get_filepath(filename))
  return dd["Y"], dd["D"], dd["M"], dd["K"], dd["Bp"], dd["Bf"], dd["A"]
end

include("RCLLadders.jl")

export gugercin_pH_msd_chain, poro_elasticity_model

end # module
