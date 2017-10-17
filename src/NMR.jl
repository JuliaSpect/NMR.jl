# Starting point for NMR work in Julia

module NMR

using Interpolations

struct ProcessedSpectrum
    re_ft :: Vector{Int32}
    im_ft :: Vector{Int32}
    params :: Dict{AbstractString, Any}
    intrng :: Vector{Tuple{Float64, Float64}}
end

Base.getindex(p::ProcessedSpectrum, param::AbstractString) = p.params[param]
Base.getindex(p::ProcessedSpectrum, ::Colon) = p.re_ft
Base.getindex(p::ProcessedSpectrum, r::Range) = p.re_ft[r]

struct Spectrum
    fid :: Vector{Int32}
    acqu :: Dict{AbstractString, Any}
    procs :: Dict{Int,ProcessedSpectrum}
    default_proc :: Int
end

Base.getindex(s::Spectrum, i::Int) = s.procs[i]
Base.getindex(s::Spectrum, ::Colon) = s[s.default_proc].re_ft
Base.getindex(s::Spectrum, r::Range) = s[s.default_proc].re_ft[r]
function Base.getindex(s::Spectrum, param::AbstractString)
    try
        s.acqu[param]
    catch err
        s[s.default_proc][param]
    end
end

Spectrum(fid :: Vector{Int32}, acqu :: Dict{Any, Any}, proc :: ProcessedSpectrum) = Spectrum(fid, acqu, Dict(1=>proc), default_proc = 1)

include("analysis.jl")
include("bruker.jl")
include("interpolation.jl")
include("utils.jl")
include("intrng.jl")


export read_bruker_binary, acq_params, interpolate_spect, analyze_lsq

end
