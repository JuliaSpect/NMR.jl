# Bruker I/O functions

"""
    read_bruker_binary(fname)

Bruker FID, and processed frequency-domain data files (fid, 1r, 1i) contain a flat
binary representation of their data points, each of which is an Int32.
"""
function read_bruker_binary(fname)
    reinterpret(Int32, open(read, fname))
end

function parse_float_list(m)
    lines = split(m)
    # The first line of a match is "(0..n)"; we don't need this.
    [parse(Float64, s) for s in split(m)[2:end]]
end

function parse_or(::Type{T}, s, default::T) where T
    v = tryparse(T, s)
    v === nothing ? default : v
end

filters = [ ( Set(["SW", "O1", "SFO1", "SF", "BF1"]),
              s -> parse(Float64, s) ),
            ( Set(["TD", "NS", "DS", "SI", "DATE"]),
              s -> parse(Int, s) ),
            ( Set(["D", "P"]),
              parse_float_list),
            ( Set(["PULPROG", "EXP"]),
              x -> strip(x)[2:end-1] ),
            ( Set(["AUTOPOS"]),
              x -> parse_or(Int, strip(x)[2:end-1], 0) )
]

function read_params(file)
    contents = read(file, String)
    matches = eachmatch(r"##\$?(.*?)=\s?([^#]+)"s, contents)
    res = Dict(m.captures[1] => parse_param(m.captures[1], m.captures[2]) for m in matches)

    # a few little tweaks
    if "O1" in keys(res) && !("O1P" in keys(res))
        res["O1P"] = res["O1"] / res["SFO1"]
    end
    res
end

function read_intrng(file)
    lines = try
        readlines(file)
    catch
        return []
    end
    if length(lines) > 2 && strip(lines[1])[1] == 'A'
        [tuple(map(s -> parse(Float64, s), split(line)[1:2])...) for line in lines[3:end]]
    elseif length(lines) > 1 && strip(lines[1])[1] == 'P'
        [tuple(map(s -> parse(Float64, s), split(line))...) for line in lines[2:end]]
    else
        []
    end
end

function parse_param(param, val)
    for (names, fun) in filters
        if param in names
            return fun(val)
        end
    end
    return strip(val)
end

"""
    ProcessedSpectrum(path :: AbstractString, procno :: Int)

Construct ProcessedSpectrum object parsing data in `procno` within `path`.
"""
function ProcessedSpectrum(path :: AbstractString, procno :: Int)
    re_ft = float(read_bruker_binary(joinpath(path, "1r")))
    im_ft = float(read_bruker_binary(joinpath(path, "1i")))
    params = read_params(joinpath(path, "procs"))
    title = read(joinpath(path, "title"), String)
    intrng = read_intrng(joinpath(path, "intrng"))
    ProcessedSpectrum(re_ft, im_ft, params, intrng, procno, title)
end

"""
    ProcessedSpectrum(path::AbstractString)

If no procno given, it is extracted from basename of `path`.
"""
ProcessedSpectrum(path::AbstractString) = ProcessedSpectrum(path, parse(Int, basename(path)))

"""
    Spectrum(path :: AbstractString, procnos :: AbstractArray{Int}, default_proc :: Int)

Use data only in the array of `procnos` within `path`, with defined `default_proc`.
"""
Spectrum(path :: AbstractString, procnos :: AbstractArray{Int}, default_proc :: Int) = begin
    fid = float(read_bruker_binary(joinpath(path, "fid")))
    acqu = read_params(joinpath(path, "acqus"))
    name = basename(dirname(path))
    expno = parse(Int, basename(path))
    procs = Dict()
    for procno in procnos
        proc_path = joinpath(path, "pdata", string(procno))
        procs[procno] = ProcessedSpectrum(proc_path, procno)
    end
    Spectrum(fid, acqu, procs, default_proc, name, expno)
end

"""
    Spectrum(path :: AbstractString, procno :: Int)

Use data only from `procno` within `path`.
"""
Spectrum(path :: AbstractString, procno :: Int) = Spectrum(path, [procno], procno)

"""
    Spectrum(path :: AbstractString, procnos :: AbstractArray{Int})

Use data only in the array of `procnos` within `path`, setting the lowest in `procnos` as the default.
"""
Spectrum(path :: AbstractString, procnos :: AbstractArray{Int}) = Spectrum(path, procnos, minimum(procnos))

"""
    Spectrum(path :: AbstractString)

Use data from every processed data file within `path`.
"""
Spectrum(path :: AbstractString) = begin
    procnos = [parse(Int,n) for n in readdir(joinpath(path, "pdata"))]
    Spectrum(path, procnos)
end
