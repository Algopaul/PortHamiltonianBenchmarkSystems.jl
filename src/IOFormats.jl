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

function loadHDF5(filename::String)
  fid = h5open(filename)
  d = Dict()
  for key in keys(fid)
    d[key] =  m_read(fid, key, fid[key])
  end
  return d
end

function m_read(::Any, ::Any, dset::HDF5.Dataset)
  return read(dset)
end
function m_read(fid, key, dset::HDF5.Group)
  return H5SparseMatrixCSC(fid, key) |> SparseMatrixCSC
end


function writeMAT(filename::String, data::Tuple{String, AbstractMatrix}...)
  d = Dict()
  for mat in data
    d[mat[1]] = mat[2]
  end
  matwrite(filename, d)
end

function loadMAT(filename::String)
  return matread(filename)
end

export writeHDF5, writeMAT
