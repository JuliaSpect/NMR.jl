# Continuous interpolation of discrete frequency-domain data
using Interpolations
import Interpolations: interpolate

function interpolate(s::Spectrum)
    l,h = limits(s)
    fn = extrapolate(interpolate(s[:], BSpline(Cubic(Natural())), OnGrid()), Flat())
    # Interpolate.jl only supports increasing ranges, so will have to
    # use flip `h` and `l` here and invert the axis in the next line.
    scaled = scale(fn, range(l; stop=h, length=length(s)))
    δ -> scaled[l + (h - δ)]
end

function resample(s::Spectrum, shifts::AbstractArray)
    intp = interpolate(s)
    intp.(shifts)
end

