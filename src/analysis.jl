"""     analyze_lsq(mix::Vector{T}, library::Vector{Vector{S}})
Performs least squares fitting of the spectrum of mixture to a library of
known spectra. Sensitive to chemical shift drifts - all spectra need to be
referenced."""
function analyze_lsq(mix::Vector{T}, library::Vector{Vector{S}}) where T<:Number where S<:Number
    lib_matrix =hcat(library...)
    guess =  lib_matrix \ mix
    residue = norm(lib_matrix * guess - mix, 2) / norm(mix, 2)
    (guess, residue)
end

function analyze_tsq(mix::Vector{T}, library::Vector{Vector{S}}) where T<:Number where S<:Number
    lib_matrix =hcat(library...)
    guess =  lib_matrix \ mix
    residue = norm(lib_matrix * guess - mix, 2) / norm(mix, 2)
    (guess, residue)
end

function analyze_lsq(mix::Spectrum, library::AbstractArray{Spectrum, 1}, dark_areas::AbstractArray{UnitRange, 1} = [])
    lib_reals = (s.procs[1].re_ft for s in library)
    lib_matrix = hcat(library...)
    guess =  lib_matrix \ mix
    residue = norm(lib_matrix * guess - mix, 2) / norm(mix, 2)
    (guess, residue)
end

