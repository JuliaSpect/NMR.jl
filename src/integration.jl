function intrng_indices(s::Spectrum)
    intrng = s[s.default_proc].intrng
    [ppmtoindex(s,i[1]):ppmtoindex(s,i[2]) for i in intrng]
end

function intrng_data(s::Spectrum)
    [s[r] for r in intrng_indices(s)]
end

function remove_rng(s::Spectrum, δ::Float64)
    remove_rng(s[s.default_proc], δ)
end

function remove_rng(p::ProcessedSpectrum, δ)
    for (i,(hi,lo)) in enumerate(p.intrng)
        if lo < δ < hi
            deleteat!(p.intrng, i)
        end
    end
end

integrate(p::ProcessedSpectrum, r::Range) = sum(p[r])

function integrate(p::ProcessedSpectrum, ref_rng::Int)
    rngs = p.intrng
    ref_int = integrate(p, rngs[ref_rng])
    [integrate(p, r)/ref_int for r in p.intrng]
end

function integrate(p::ProcessedSpectrum)
    [integrate(p, r)/ref_int for r in p.intrng]
end

integrate(s::Spectrum, args...) = integrate(s[s.default_proc], args...)