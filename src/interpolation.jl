# Continuous interpolation of discrete frequency-domain data
using Interpolations
import Interpolations: interpolate

"""
    interpolate(s::Spectrum)

Create an interpolate objection with data and limits of Spectrum `s`
"""
function interpolate(s::Spectrum)
    l,h = limits(s)
    fn = extrapolate(interpolate(s[:], BSpline(Cubic(Natural(OnGrid())))), Flat())
    scaled = scale(fn, range(l; stop=h, length=length(s)))
    δ -> scaled(h + (l - δ))
end

"""
    resample(s::Spectrum, shifts::AbstractArray)

Interpolate data from Spectrum `s` and resample data points at
chemical shifts defined by `shifts`.
"""
function resample(s::Spectrum, shifts::AbstractArray)
    intp = interpolate(s)
    intp.(shifts)
end

