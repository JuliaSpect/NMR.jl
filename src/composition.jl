#using Interpolations

import Base: +

function +(s1::Spectrum, s2::Spectrum)
    l1,h1 = limits(s1)
    l2,h2 = limits(s2)
    res = min(freq_resolution(s1), freq_resolution(s2))
    l = min(l1, l2)
    h = max(h1, h2)
    shifts = l:res:h
    vals = similar(shifts, Float64)
    re1, re2 = interpolate(s1), interpolate(s2)
    for (i,δ) in enumerate(shifts)
        vals[i] = re1[δ] + re2[δ]
    end
    # TODO: make a proper processed Spectrum
    (shifts, vals)
end
