#using Interpolations

import Base: +


function +(s1::Spectrum, s2::Spectrum)
    shifts = union_shifts([s1, s2])
    vals = resample(s1, shifts) .+ resample(s2, shifts)
    # TODO: make a proper processed Spectrum
    (shifts, vals)
end
