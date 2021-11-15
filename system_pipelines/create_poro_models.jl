using PortHamiltonianBenchmarkSystems, MAT

for n in [320, 980, 1805]
  E, J, R, B = poro_elasticity_model(n=n)
  filename = "PoroModelN$(2n_cells)"
  writeHDF5(filename*".h5"; E, J, R, B)
  writeMAT(filename*".mat"; E, J, R, B)
end
