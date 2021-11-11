using HDF5, H5Sparse, MAT

function writeHDF5(filename::String, data::Tuple{String, AbstractMatrix}...)
  fid = h5open(filename, "cw")
  for mat in data
    m_write_hdf5(filename, mat[1], mat[2])
  end
  close(fid)
end

m_write_hdf5(f::String, varname::String, var::Matrix) = h5write(f, varname, var)
m_write_hdf5(f::String, varname::String, var::SparseMatrixCSC) = H5SparseMatrixCSC(f, varname, var)

function writeMAT(filename::String, data::Tuple{String, AbstractMatrix}...)
  d = Dict()
  for mat in data
    d[mat[1]] = mat[2]
  end
  matwrite(filename, d)
end

export writeHDF5, writeMAT
