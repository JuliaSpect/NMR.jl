limits(s::Spectrum) = limits(s["O1P"], s["SW"], s["SF"], s["BF1"])
limits(o1p, sw, sf, bf) = begin
    shift = 1e6(bf-sf)/bf
    (o1p - sw/2 + shift, o1p + sw/2 + shift)
end

function chemical_shifts(s::Spectrum)
    lo,hi = limits(s)
    linspace(hi,lo,length(s[:]))
end


"""     toppm(f, bf)
Hz to ppm conversion of a function
f: Function that takes frequency. (Hz)
bf: Base frequency. (Hz)
δ: Chemical shift. (ppm)"""
function toppm(f, bf)
    δ -> f(bf*(1.0+δ*1e-6))
end

function sr(sf, bf)
    return 1e6(sf-bf)
end

function ppmtoindex(δ, sw, o1p, sr, sf)
    max_shift = o1p + 0.5sw - sr/sf 
end

function ppmtoindex(s::Spectrum, δ)
    min_δ ,max_δ = limits(s)
    Int(round((max_δ - δ)/(max_δ - min_δ) * s["SI"]))
end

function ppmtoindex(s::Spectrum, rng::Tuple{Float64,Float64})
    r1,r2 = rng
    if r1>r2
        ppmtoindex(s,r1):ppmtoindex(s,r2)
    else
        ppmtoindex(s,r2):ppmtoindex(s,r1)
    end
end

Base.length(s::Spectrum) = length(s[:])
Base.length(p::ProcessedSpectrum) = length(p[:])
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
