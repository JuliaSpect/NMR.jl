"""
    intrng(s::Spectrum)

Return the integral ranges defined in the default processed spectrum in `s`.
"""
intrng(s::Spectrum) = s[s.default_proc].intrng

"""
    intrng_indices(s::Spectrum)

Return the integral ranges defined in `s` as an array of ranges of indices.
"""
function intrng_indices(s::Spectrum)
    rng = intrng(s)
    [ppmtoindex(s,i[1]):ppmtoindex(s,i[2]) for i in rng]
end

"""
    intrng_data(s::Spectrum)

Return the data encompassed by integral ranges defined in `s` as separate arrays.
"""
function intrng_data(s::Spectrum)
    [@view s[r] for r in intrng_indices(s)]
end

"""
    intrng_shifts(s::Spectrum)

Return the integral ranges defined in `s` as an array of ranges in ppm.
"""
function intrng_shifts(s::Spectrum)
    rng = intrng(s)
    [range(i[1]; stop=i[2], length=ppmtoindex(s,i[2])-ppmtoindex(s,i[1])+1) for i in rng]
end

function find_rng(p::ProcessedSpectrum, δ::Float64)
    for (i,(hi,lo)) in enumerate(p.intrng)
        if lo < δ < hi
            return i
        end
    end
end

find_rng(s::Spectrum, δ::Float64) = find_rng(s[s.default_proc], δ)

function remove_rng(sp::Union{Spectrum,ProcessedSpectrum}, δ::Float64)
    t = deepcopy(sp)
    remove_rng!(t, δ)
    t
end

remove_rng!(s::Spectrum, δ::Float64) = remove_rng!(s[s.default_proc], δ)
remove_rng!(p::ProcessedSpectrum, δ::Float64) = remove_rng!(p, find_rng(p, δ))
remove_rng!(p::ProcessedSpectrum, n::Int) = deleteat!(p.intrng, n)
remove_rng!(p::ProcessedSpectrum, x::Nothing) = nothing

"""
    integrate(v::Vector, r::UnitRange)

Return the integral of the vector `v` within the unit range `r`.
"""
integrate(v::Vector, r::UnitRange) = sum(v[r])

"""
    integrate(p::ProcessedSpectrum, r::UnitRange)

For processed spectrum `p`.
"""
integrate(p::ProcessedSpectrum, r::UnitRange) = integrate(p[:], r)

"""
    integrate(s::Spectrum, r::UnitRange)

For NMR spectrum `s`, using the default processed spectrum as defined within `s`.
"""
integrate(s::Spectrum, r::UnitRange) = integrate(s[s.default_proc], r)

"""
    integrate(s::Spectrum, ppm_range::Tuple{Float64,Float64})

For NMR spectrum `s` within `ppm_range`.
"""
function integrate(s::Spectrum, ppm_range::Tuple{Float64,Float64})
    r1 = ppmtoindex(s, ppm_range[1])
    r2 = ppmtoindex(s, ppm_range[2])
    integrate(s, r1:r2)
end

"""
    integrate(s::Spectrum)

If no range is specified, return the integrals of every integral range defined in `s`.
"""
function integrate(s::Spectrum)
    [integrate(s, r) for r in intrng(s)]
end

"""
    integrate(s::Spectrum, ref_rng::Int)

Return the integrals of every integral range defined in `s`, normalizing to the `ref_rng`th integral.
"""
function integrate(s::Spectrum, ref_rng::Int)
    rngs = intrng(s)
    ref_int = integrate(s, rngs[ref_rng])
    integrate(s)./ref_int
end

"""
    integrate(s::Spectrum, δ::Float64)

Return the integrals of every integral range defined in 's', normalizing to the
integral range encompassing the chemical shift `δ`.
"""
integrate(s::Spectrum, δ::Float64) = integrate(s, find_rng(s, δ))
