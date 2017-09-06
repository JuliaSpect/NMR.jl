# Bruker I/O functions

"""
    read_bruker_binary(fname)
Bruker FID, and processed frequency-domain data files (fid, 1r, 1i) contain a flat
binary representation of their data points, each of which is an Int32."""
function read_bruker_binary(fname)
    reinterpret(Int32, open(read, fname))
end

filters = [ ( Set(["SW", "O1", "SFO1"]),
              float )
          ]

function acq_params(acqufile)
    matches = eachmatch(r"##\$(.*)=\s?(.*)"m, open(readstring, acqufile))
    res = Dict(m.captures[1] => parse_param(m.captures[1], m.captures[2]) for m in matches)

    # a few little tweaks
    res["O1P"] = res["O1"] / res["SFO1"]
    res
end

function parse_param(param, val)
    for (names, fun) in filters
        if param in names
            return fun(val)
        end
    end
    return strip(val)
end
