# Continuous interpolation of discrete frequency-domain data

"""
    interpolate_spectrum(data, cs_range)"""
function interpolate_spectrum(data, cs_range, options... = [BSpline(Cubic(Natural())), OnGrid()])
    lo, hi = cs_range
    N = length(data)
    f = interpolate(data, options...)
    # lo -> N, hi -> 1
    d -> f[N - (d - lo)/(hi - lo)*(N - 1)]
end
