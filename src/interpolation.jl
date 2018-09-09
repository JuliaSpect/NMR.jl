# Continuous interpolation of discrete frequency-domain data
using Interpolations
import Interpolations: interpolate

function interpolate(s::Spectrum)
    l,h = limits(s)
    fn = extrapolate(interpolate(s[:], BSpline(Cubic(Natural())), OnGrid()), Flat())
    scaled = scale(fn, range(h; stop=l, length=length(s)))
    δ -> scaled[δ]
end

function resample(s::Spectrum, shifts::AbstractArray)
    intp = interpolate(s)
    intp.(shifts)
end

