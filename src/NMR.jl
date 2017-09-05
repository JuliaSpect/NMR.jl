# Starting point for NMR work in Julia

module NMR

using Interpolations

include("analysis.jl")
include("bruker.jl")
include("interpolation.jl")


export read_bruker_binary, acq_params, interpolate_spect, analyze_lsq

end
