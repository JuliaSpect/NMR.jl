# Starting point for NMR work in Julia

module NMR

export Spectrum, ProcessedSpectrum, plot, plot!

const Intrng = Tuple{Float64, Float64}

struct ProcessedSpectrum
    re_ft :: Vector{Float64}
    im_ft :: Vector{Float64}
    params :: Dict{String, Any}
    intrng :: Array{Intrng,1}
    procno :: Int
    title :: String
end
Base.getindex(p::ProcessedSpectrum, n::Int) = p.intrng[n]
Base.getindex(p::ProcessedSpectrum, param::AbstractString) = p.params[param]
Base.getindex(p::ProcessedSpectrum, ::Colon) = p.re_ft
Base.getindex(p::ProcessedSpectrum, r::Range) = p.re_ft[r]

mutable struct Spectrum
    fid :: Vector{Float64}
    acqu :: Dict{String, Any}
    procs :: Dict{Int,ProcessedSpectrum}
    default_proc :: Int
    name :: String
    expno :: Int
end

Base.getindex(s::Spectrum, i::Int) = s.procs[i]
Base.getindex(s::Spectrum, ::Colon) = s.procs[s.default_proc].re_ft
Base.getindex(s::Spectrum, r::Range) = s.procs[s.default_proc].re_ft[r]
Base.getindex(s::Spectrum, rng::Tuple{Float64,Float64}) = s[ppmtoindex(s,rng)]
Base.getindex(s::Spectrum, δ::Float64) = s[s.default_proc].re_ft[ppmtoindex(s, δ)]
function Base.getindex(s::Spectrum, param::AbstractString)
    try
        s.acqu[param]
    catch err
        s[s.default_proc][param]
    end
end
Base.setindex!(s::Spectrum, d::AbstractArray, ::Colon) = (s[s.default_proc].re_ft .= d)
Base.setindex!(s::Spectrum, d::AbstractArray, r::Range) = (s[s.default_proc].re_ft[r] .= d)
Base.setindex!(s::Spectrum, d, rng::Tuple{Float64, Float64}) = (s[s.default_proc].re_ft[ppmtoindex(s,rng)]=d)

Spectrum(fid :: Vector{Float64}, acqu :: Dict{Any, Any}, proc :: ProcessedSpectrum) = Spectrum(fid, acqu, Dict(1=>proc), 1, "", "")

include("baseline.jl")
include(joinpath("bruker", "read.jl"))
include(joinpath("bruker", "write.jl"))
include("composition.jl")
include("decomposition.jl")
include("integration.jl")
include("interpolation.jl")
include("plotting.jl")
include("utils.jl")

export read_bruker_binary, acq_params, interpolate_spect, lsq_analyze,
       overlay!, candidates, cand_signals, guess_matrices, limits, chemical_shifts,
       Spectrum, ProcessedSpectrum

end
