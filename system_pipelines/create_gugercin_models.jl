using PortHamiltonianBenchmarkSystems, MAT

for n_cells in [50, 500, 10_000]
  J, R, Q, B = gugercin_pH_msd_chain(n_cells = n_cells)
  filename = "GugercinN$(2n_cells)"
  writeHDF5(filename*".h5"; J, R, Q, B)
  writeMAT(filename*".mat"; J, R, Q, B)
end
