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