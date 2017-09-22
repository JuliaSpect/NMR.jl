# Starting point for NMR work in Julia

module NMR

using Interpolations

struct ProcessedSpectrum
    re_ft :: Vector{Int32}
    im_ft :: Vector{Int32}
    params :: Dict{Any, Any}
end

struct Spectrum
    fid :: Vector{Int32}
    acqu :: Dict{Any, Any}
    procs :: Dict{Int,ProcessedSpectrum}
    default_proc :: Int
end

Spectrum(fid :: Vector{Int32}, acqu :: Dict{Any, Any}, proc :: ProcessedSpectrum) = Spectrum(fid, acqu, Dict(1=>proc), default_proc = 1)

include("analysis.jl")
include("bruker.jl")
include("interpolation.jl")
include("utils.jl")


export read_bruker_binary, acq_params, interpolate_spect, analyze_lsq

end
