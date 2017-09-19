# Bruker I/O functions

"""     read_bruker_binary(fname)
Bruker FID, and processed frequency-domain data files (fid, 1r, 1i) contain a flat
binary representation of their data points, each of which is an Int32."""
function read_bruker_binary(fname)
    reinterpret(Int32, open(read, fname))
end

function parse_float_list(m)
    lines = split(m)
    # The first line of a match is "(0..n)"; we don't need this.
    [float(s) for s in split(m)[2:end]]
end

filters = [ ( Set(["SW", "O1", "SFO1"]),
              float ),
            ( Set(["TD", "NS", "DS", "SI"]),
              s -> parse(Int, s) ),
            ( Set(["D", "P"]),
              parse_float_list),
            ( Set(["PULPROG"]),
              x -> strip(x)[2:end-1])
]

function read_params(file)
    matches = eachmatch(r"##\$?(.*?)=\s?([^#]+)"s, open(readstring, file))
    res = Dict(m.captures[1] => parse_param(m.captures[1], m.captures[2]) for m in matches)

    # a few little tweaks
    if "O1" in keys(res) && !("O1P" in keys(res))
        res["O1P"] = res["O1"] / res["SFO1"]
    end
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

Spectrum(path :: AbstractString) = begin
    fid = read_bruker_binary(joinpath(path, "fid"))
    acqu = read_params(joinpath(path, "acqu"))
    procs = Dict()
    for procno in readdir(joinpath(path, "pdata"))
        proc_path = joinpath(path, "pdata", procno)
        re_ft = read_bruker_binary(joinpath(proc_path, "1r"))
        im_ft = read_bruker_binary(joinpath(proc_path, "1i"))
        params = read_params(joinpath(proc_path, "proc"))
        procs[parse(Int, procno)] = ProcessedSpectrum(re_ft, im_ft, params)
    end
    Spectrum(fid, acqu, procs)
end

Spectrum(path :: AbstractString, procnos :: AbstractArray{T} where T <: Integer) = begin
    fid = read_bruker_binary(joinpath(path, "fid"))
    acqu = read_params(joinpath(path, "acqu"))
    procs = Dict()
    for procno in procnos
        proc_path = joinpath(path, "pdata", string(procno))
        re_ft = read_bruker_binary(joinpath(proc_path, "1r"))
        im_ft = read_bruker_binary(joinpath(proc_path, "1i"))
        params = read_params(joinpath(proc_path, "proc"))
        procs[procno] = ProcessedSpectrum(re_ft, im_ft, params)
    end
    Spectrum(fid, acqu, procs)
end
