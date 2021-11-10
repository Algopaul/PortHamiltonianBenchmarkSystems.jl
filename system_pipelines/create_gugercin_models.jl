using PortHamiltonianBenchmarkSystems, HDF5

for n_cells in [50, 500, 10_000]
  J, R, Q, B = gugercin_pH_msd_chain(n_cells = n_cells)
  h5open("GugercinN$(2*n_cells).h5", "w") do file
    write(file, "J", J, "R", Q, "Q", B, "B")
  end
end
