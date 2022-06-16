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

export poro_elasticity_model
