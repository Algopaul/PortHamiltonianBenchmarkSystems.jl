using MD5

hex(s) = string(s, base=16)
"""
Downloads a dataset located at a given url to a path. The function is taken
from Flux.jl/Data.jl
``https://github.com/FluxML/Flux.jl/blob/ea26f45a1f4e93d91b1e8942c807f8bf229d5775/src/data/Data.jl#L26``
"""
function download_dataset(url, path, hash=nothing)
  println("Downloading...")
  if hash !== nothing
    tmppath = tempname()
    download(url, tmppath)
    hash_download = bytes2hex(open(md5, tmppath))
    if hash_download !== hash
      msg  = "Hash Mismatch!\n"
      msg *= "  Expected md5: $hash\n"
      msg *= "  Calculated md5: $hash_download\n"
      error(msg)
    end
    mv(tmppath, path; force = true)
  else
    download(url, path)
  end
end

function is_new(path::String, online_hash)
  if online_hash === nothing
    @warn "No hash-value provided. File might have been modified. You can force_download trigger a new download."
    return false
  else
    offline_hash = bytes2hex(open(md5, path))
    return offline_hash !== online_hash
  end
end

function get_filepath(filename)
  path = @__DIR__
  path *= "/../data/"
  return path*"$filename"
end

function download_system_data_if_required(filename, url, hash=nothing; force_download = false)
  filepath = get_filepath(filename)
  if !isfile(filepath)
    download_dataset(url, filepath, hash)
  elseif is_new(filepath, hash)
    println("Local file is modified.")
    download_dataset(url, filepath, hash)
  elseif force_download === true
    download_dataset(url, filepath, hash)
  end
end
