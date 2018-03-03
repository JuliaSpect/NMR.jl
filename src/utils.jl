limits(s::Spectrum) = limits(s["O1P"], s["SW"], s["SF"], s["BF1"])
limits(o1p, sw, sf, bf) = begin
    shift = 1e6(bf-sf)/bf
    (o1p - sw/2 + shift, o1p + sw/2 + shift)
end

function within(s::Spectrum, δ)
    l,h = limits(s)
    l < δ < h
end

function chemical_shifts(s::Spectrum)
    lo,hi = limits(s)
    linspace(hi,lo,length(s[:]))
end


"""     hztoppm(f, bf)
Hz to ppm conversion of a function
f: Function that takes frequency. (Hz)
bf: Base frequency. (Hz)
δ: Chemical shift. (ppm)"""
function hztoppm(f, bf)
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
    @. Int(cld(s["SI"]*(max_δ - δ), (max_δ - min_δ)))
end

function ppmtoindex(s::Spectrum, rng::Tuple{Float64,Float64})
    r1,r2 = rng
    if r1>r2
        ppmtoindex(s,r1):ppmtoindex(s,r2)
    else
        ppmtoindex(s,r2):ppmtoindex(s,r1)
    end
end

hztoindex(f, sw, sf, si) = Int(cld(f*si, sw*sf))
hztoindex(s::Spectrum, f) = hztoindex(f, s["SW"], s["SF"], s["SI"])

Base.length(s::Spectrum) = length(s[:])
Base.length(p::ProcessedSpectrum) = length(p[:])

function union_range(ss::AbstractArray{Spectrum})
    lims = [limits(s) for s in ss]
    l = minimum(lim[1] for lim in lims)
    h = maximum(lim[2] for lim in lims)
    h,l
end

function union_shifts(ss::AbstractArray{Spectrum})
    h,l = union_range(ss)
    res = minimum(freq_resolution(s) for s in ss)
    h:-res:l
end

union_range(s::Spectrum) = union_range([s])
union_shifts(s::Spectrum) = union_shifts([s])

freq_resolution(s::Spectrum) = s["SW"] / length(s)

title(s::Spectrum) = s[s.default_proc].title

# Returns a copy of s with all but the given
# ranges zeroed out. The default proc for s will
# have its intrng adjusted to Δs as well.
function extract(s::Spectrum, Δs::AbstractArray{Intrng})
    res = deepcopy(s)
    res[res.default_proc].re_ft = zeros(length(s))
    res[res.default_proc].im_ft = zeros(length(s))
    res[res.default_proc].intrng = Δs
    for Δ in Δs
        rng = ppmtoindex(s, Δ)
        res[res.default_proc].re_ft[rng] .= s[s.default_proc].re_ft[rng]
        res[res.default_proc].im_ft[rng] .= s[s.default_proc].im_ft[rng]
    end
    res
end

extract(s::Spectrum, Δ::Intrng) = extract(s, [Δ])

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
