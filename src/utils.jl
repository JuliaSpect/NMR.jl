limits(s :: Spectrum) = limits(s.acqu["O1P"], s.acqu["SW"])
limits(o1p, sw) = (o1p - sw/2, o1p + sw/2)

"""     toppm(f, bf)
Hz to ppm conversion of a function
f: Function that takes frequency. (Hz)
bf: Base frequency. (Hz)
δ: Chemical shift. (ppm)"""
function toppm(f, bf)
    δ -> f(bf*(1.0+δ*1e-6))
end

### Plotting functions

import Plots: plot

function plot(s :: Spectrum, fid :: Bool = false)
    if !fid
    lo,hi = limits(s)
    shifts = linspace(lo, hi, length(s.re_ft))    
    # (shifts, s.re_ft)
    plot(shifts, s.re_ft, xflip=true)
    end
end

### Pulse power profile

"""     powerprofile(pulse, dt, sfo1)
Pulse power profile for arbitrary pulse shape.
pulse: Amplitude profile of pulse.
dt: 'Sampling frequency' — unit of time between points in time. (s)
sfo1: Centre frequency of pulse. (Hz)"""
function powerprofile(pulse, dt, sfo1)
    ft = abs2.(fft(pulse))
    l = length(pulse)
    df = 1/(l*dt)
    fn = interpolate(ft[1:div(l,2)], BSpline(Cubic(Natural())), OnGrid())
    max = fn[1]
    f -> fn[abs(f-sfo1)/df+1]/max
end
