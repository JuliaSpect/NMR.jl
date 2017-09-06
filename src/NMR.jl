# Starting point for NMR work in Julia

module NMR

using Interpolations

struct Spectrum
    fid :: Vector{Int32}
    re_ft :: Vector{Int32}
    im_ft :: Vector{Int32}
    acqu :: Dict{Any, Any}
end

include("analysis.jl")
include("bruker.jl")
include("interpolation.jl")
include("utils.jl")


export read_bruker_binary, acq_params, interpolate_spect, analyze_lsq

end
