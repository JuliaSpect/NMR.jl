__precompile__(true)

module NMR

const Intrng = Tuple{Float64, Float64}

mutable struct ProcessedSpectrum
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
Base.getindex(p::ProcessedSpectrum, a::AbstractArray) = p.re_ft[a]

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
Base.getindex(s::Spectrum, a::AbstractArray) = s.procs[s.default_proc].re_ft[a]
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
include("bruker/read.jl")
include("bruker/write.jl")
include("composition.jl")
include("decomposition.jl")
include("integration.jl")
include("interpolation.jl")
include("plotting.jl")
include("utils.jl")

export Spectrum, ProcessedSpectrum, dump, # constructors, I/O
       plot, plot!, # plotting
       lsq_analyze, candidates, decompose, # decomposition
       ppmtoindex, hztoindex, hztoppm, # unit conversion
       baseline_correct!, # processing
       limits, chemical_shifts, extract, copy! # utility functions
       interpolate, resample # interpolation
end
