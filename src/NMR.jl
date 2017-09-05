module NMR
export read_bruker_binary
using Interpolations

# Starting point for NMR work in Julia

# Bruker I/O functions
"""
    read_bruker_binary(fname)
Bruker FID, and processed frequency-domain data files (fid, 1r, 1i) contain a flat
binary representation of their data points, each of which is an Int32."""
function read_bruker_binary(fname)
    reinterpret(Int32, open(read, fname))
end

"""
    interpolate_spect(data, cs_range)"""
function interpolate_spect(data, cs_range)
    lo, hi = cs_range
    N = length(data)
    f = interpolate(data, BSpline(Cubic(Natural())), OnGrid())
    # lo -> N, hi -> 1
    d -> f[N - (d - lo)/(hi - lo)*(N - 1)]
end

function acq_params(acqufile)
    matches = eachmatch(r"##\$(.*)=\s?(.*)", open(readstring, acqufile))
    
end

end
