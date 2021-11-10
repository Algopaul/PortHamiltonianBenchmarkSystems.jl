using PortHamiltonianBenchmarkSystems, MAT

for n_cells in [50, 500, 10_000]
  J, R, Q, B = PortHamiltonianBenchmarkSystems.gugercin_pH_msd_chain(n_cells = n_cells)
  matwrite("GugercinN$(2*n_cells).mat", Dict(
       "J" => J,
       "R" => R,
       "Q" => Q,
       "B" => B,
      ))
end
