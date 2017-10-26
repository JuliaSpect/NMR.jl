# Starting point for NMR work in Julia

module NMR

export Spectrum, ProcessedSpectrum, plot, plot!

struct ProcessedSpectrum
    re_ft :: Vector
    im_ft :: Vector
    params :: Dict{AbstractString, Any}
    intrng :: Vector{Tuple{Float64, Float64}}
    procno :: Int
end

Base.getindex(p::ProcessedSpectrum, param::AbstractString) = p.params[param]
Base.getindex(p::ProcessedSpectrum, ::Colon) = p.re_ft
Base.getindex(p::ProcessedSpectrum, r::Range) = p.re_ft[r]

mutable struct Spectrum
    fid :: Vector{Int32}
    acqu :: Dict{AbstractString, Any}
    procs :: Dict{Int,ProcessedSpectrum}
    default_proc :: Int
    name :: AbstractString
    expno :: Int
end

Base.getindex(s::Spectrum, i::Int) = s.procs[i]
Base.getindex(s::Spectrum, ::Colon) = s[s.default_proc].re_ft
Base.getindex(s::Spectrum, r::Range) = s[s.default_proc].re_ft[r]
Base.getindex(s::Spectrum, rng::Tuple{Float64,Float64}) = s[ppmtoindex(s,rng)]
function Base.getindex(s::Spectrum, param::AbstractString)
    try
        s.acqu[param]
    catch err
        s[s.default_proc][param]
    end
end
Base.setindex!(s::Spectrum, d::AbstractArray, ::Colon) = (s[s.default_proc].re_ft .= d)

Spectrum(fid :: Vector{Int32}, acqu :: Dict{Any, Any}, proc :: ProcessedSpectrum) = Spectrum(fid, acqu, Dict(1=>proc), 1, "", "")

include("baseline.jl")
include("bruker.jl")
include("composition.jl")
include("decomposition.jl")
include("integration.jl")
include("interpolation.jl")
include("plotting.jl")
include("utils.jl")

export read_bruker_binary, acq_params, interpolate_spect, lsq_analyze,
       overlay!, candidates, cand_signals, guess_matrices, limits, chemical_shifts

end
