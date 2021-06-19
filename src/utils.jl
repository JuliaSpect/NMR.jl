# Utility functions, mostly for unit conversions.
# NMR quantities in Bruker notation [1]:
# BF1: Basic transmitter frequency (MHz)
# SF: Spectrometer frequency (MHz)
# SR: Spectrum reference frequency (Hz)
# SR = (SF - BF1) * 10^6
# SFO1: Transmitter frequency (MHz)
# O1: Transmitter frequency offset (Hz)
# O1 = (SFO1 - BF1) * 10^6
# Ω: Absolute frequency (MHz)
# ω: Frequency relative to SF1 (Hz)
# δ: Chemical shift δ (ppm)
# δ = 10^6 * (Ω - SF) / SF = ω / SF
# 1. Numerical indices indicate different nuclei

import Base: in

"""
    limits(s::Spectrum)

Return the chemical shift limits of the Spectrum `s` in units of ppm.
"""
limits(s::Spectrum) = limits(s["O1P"], s["SW"], s["SF"], s["BF1"])

"""
    limits(o1p, sw, sf, bf)

Return the chemical shift limits given `o1p`, `sw`, `sf`, `bf`.
"""
limits(o1p, sw, sf, bf) = begin
    shift = 1e6(bf-sf)/bf
    (o1p - sw/2 + shift, o1p + sw/2 + shift)
end

function in(δ, s::Spectrum)
    l,h = limits(s)
    l < δ < h
end


"""
    chemical_shifts(s::Spectrum)

Return an iterator counting the chemical shift in ppm for each data point in Spectrum `s`
"""
function chemical_shifts(s::Spectrum)
    lo,hi = limits(s)
    range(hi; stop=lo, length=length(s[:]))
end


"""
    ppmtomhz_abs(δ, bf, sr = 0.0)

Convert chemical shift `δ` to absolute frequency Ω in MHz.
`bf` in MHz and `sr` in Hz.
"""
function ppmtomhz_abs(δ, bf, sr = 0.0)
    bf = bf * 1e6
    1e-6(bf+sr)δ + bf + sr
end

"""
    ppmtohz(δ, bf, sr = 0.0)

Convert chemical shift `δ` in ppm to relative frequency ω in Hz.
`b` in MHz and `sr` in Hz.
"""
function ppmtohz(δ, bf, sr = 0.0)
    (bf + sr*1e-6)δ
end

"""
    hztoppm(ω, bf, sr = 0.0)

Convert relative frequency `ω` in Hz to chemical shift δ in ppm .
`bf` in MHz and `sr` in Hz.
"""
function hztoppm(ω, bf, sr = 0.0)
    ω / (bf + sr*1e-6)
end

"""
    sr(sf, bf)

Return spectral reference (SR in Bruker notation) in Hz.
`sf` and `bf` in MHz.
"""
sr(sf, bf) = 1e6(sf - bf)

"""
    ppmtoindex(s::Spectrum, δ)

Return the index in the spectrum `s` corresponding to
chemical shift `δ`.
"""
function ppmtoindex(s::Spectrum, δ)
    min_δ ,max_δ = limits(s)
    @. Int(cld(s["SI"]*(max_δ - δ), (max_δ - min_δ)))
end


"""
    ppmtoindex(s::Spectrum, rng::Tuple{Float64,Float64})

Return the index in the spectrum `s` corresponding to the range `rng`.
"""
function ppmtoindex(s::Spectrum, rng::Tuple{Float64,Float64})
    r1,r2 = rng
    if r1>r2
        ppmtoindex(s,r1):ppmtoindex(s,r2)
    else
        ppmtoindex(s,r2):ppmtoindex(s,r1)
    end
end

"""
    hztoindex(f, sw, sf, si)

Return the index of a spectrum with spectra width `sw`, spectrometer frequency `sf`,
number of points `si`, given the frequency `f`.
"""
hztoindex(f, sw, sf, si) = Int(cld(f*si, sw*sf))

"""
    hztoindex(s::Spectrum, f)

Return the index of a spectrum `s` given the frequency `f`.
"""
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

"""
    extract(s::Spectrum, Δs::AbstractArray{Intrng})

Return a copy of Spectrum `s` with all but the given ranges `Δs` zeroed out.
The default proc for the extracted spectrum will have its intrng adjusted to `Δs` as well.
"""
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

"""
    extract(s::Spectrum, Δ::Intrng)

For a single integral range.
"""
extract(s::Spectrum, Δ::Intrng) = extract(s, [Δ])

"""Tunes `param` until `expr` evaluates to zero within δ.
`expr` must be monotonically increasing in terms of `param`."""
macro binary_opt(expr, param, min, max, δ, nattempts=5)
    :(m = $(esc(min));M = $(esc(max));d=$(esc(δ));z=zero(d);n=$(esc(nattempts));
    for i=1:n
        $(esc(param)) = (m+M)/2
        e = $(esc(expr))
        # println("val: $e, param:$((m+M)/2), min: $m, max: $M")
        if norm(e) < d
            break
        elseif e < z
            m = $(esc(param))
        else
            M = $(esc(param))
        end
    end)
end

### Pulse power profile

"""     powerprofile(pulse, dt, sfo1)

Pulse power profile for arbitrary pulse shape.
pulse: Amplitude profile of pulse.
dt: 'Sampling frequency' — unit of time between points in time. (s)
sfo1: Centre frequency of pulse. (Hz)
"""
function powerprofile(pulse, dt, sfo1)
    ft = abs2.(fft(pulse))
    l = length(pulse)
    df = 1/(l*dt)
    fn = interpolate(ft[1:div(l,2)], BSpline(Cubic(Natural())), OnGrid())
    max = fn[1]
    f -> fn[abs(f-sfo1)/df+1]/max
end
