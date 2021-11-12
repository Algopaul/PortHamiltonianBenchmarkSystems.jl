using SHA
"""
Downloads a dataset located at a given url to a path. The function is taken
from Flux.jl/Data.jl
``https://github.com/FluxML/Flux.jl/blob/ea26f45a1f4e93d91b1e8942c807f8bf229d5775/src/data/Data.jl#L26``
"""
function download_dataset(url, path, hash=nothing)
  if hash !== nothing
    tmppath = tempname()
    download(url, tmppath)
    hash_download = open(tmppath) do f
      bytes2hex(sha256(f))
    end
    if hash_download !== hash
      msg  = "Hash Mismatch!\n"
      msg *= "  Expected sha256: $hash\n"
      msg *= "  Calculated sha256: $hash_download\n"
      error(msg)
    end
    mv(tmppath, path; force = true)
  else
    download(url, path)
  end
end
