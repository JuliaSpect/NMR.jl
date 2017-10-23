function intrng_indices(s::Spectrum)
    intrng = s[s.default_proc].intrng
    [ppmtoindex(s,i[1]):ppmtoindex(s,i[2]) for i in intrng]
end

function intrng_data(s::Spectrum)
    [s[r] for r in intrng_indices(s)]
end

function intrng_shifts(s::Spectrum) 
    intrng = s[s.default_proc].intrng
    [linspace(i[1],i[2],ppmtoindex(s,i[2])-ppmtoindex(s,i[1])+1) for i in intrng]
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

integrate(v::Vector, r::Range) = sum(v[r])

integrate(p::ProcessedSpectrum, r::Range) = integrate(p[:], r)

integrate(s::Spectrum, r::Range) = integrate(s[s.default_proc], r)

function integrate(s::Spectrum, ppm_range::Tuple{Float64,Float64})
    r1 = ppmtoindex(s, ppm_range[1])
    r2 = ppmtoindex(s, ppm_range[2])
    integrate(s, r1:r2)
end

function integrate(s::Spectrum)
    [integrate(s, r) for r in s[s.default_proc].intrng]
end

function integrate(s::Spectrum, ref_rng::Int)
    rngs = s[s.default_proc].intrng
    ref_int = integrate(s, rngs[ref_rng])
    integrate(s)./ref_int
end