# Bruker I/O functions

"""
    read_bruker_binary(fname)
Bruker FID, and processed frequency-domain data files (fid, 1r, 1i) contain a flat
binary representation of their data points, each of which is an Int32."""
function read_bruker_binary(fname)
    reinterpret(Int32, open(read, fname))
end

function acq_params(acqufile)
    matches = eachmatch(r"##\$(.*)=\s?(.*)", open(readstring, acqufile))
    
end
