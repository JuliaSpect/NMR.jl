function intrng_indices(s::Spectrum)
    intrng = s[s.default_proc].intrng
    [ppmtoindex(s,i[1]):ppmtoindex(s,i[2]) for i in intrng]
end

function intrng_data(s::Spectrum)
    [s[r] for r in intrng_indices(s)]
end
