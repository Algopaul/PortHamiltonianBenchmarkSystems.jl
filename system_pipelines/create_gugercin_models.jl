using PortHamiltonianBenchmarkSystems, MAT

for n_cells in [50, 500, 10_000]
  J, R, Q, B = PortHamiltonianBenchmarkSystems.gugercin_pH_msd_chain(n_cells = n_cells)
  filename = "GugercinN$(2n_cells)"
  writeMAT(filename*".mat", ("J", J), ("R", R), ("Q", Q), ("B", B))
  writeHDF5(filename*".h5", ("J", J), ("R", R), ("Q", Q), ("B", B))
end
