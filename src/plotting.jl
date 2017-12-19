import Plots: plot, plot!, @colorant_str, annotate!, text, font
using Formatting

function plot_limits(s::NMR.Spectrum, margin = 0.05)
    m, M = (minimum(minimum(d) for d in intrng_data(s)),
            maximum(maximum(d) for d in intrng_data(s)))
    Δ = M - m
    (m - margin*Δ, M + margin * Δ)
end

function plot_limits(a::Array{Float64,1}, margin = 0.05)
    m, M = extrema(a)
    Δ = M - m
    (m - margin*Δ, M + margin * Δ)
end

plot(s::Spectrum; integrate = false, kw...) = begin
    lo,hi = limits(s)
    shifts = linspace(hi, lo, length(s))
    p = plot(shifts, s[:], framestyle=:box)
    if integrate == true
        integration_plot(s)
    end
    plot!(xflip=true, yticks=[], leg=false; kw...)
end

function plot(spectra::AbstractArray{Spectrum}; kw...)
    shifts = union_shifts(spectra)
    labels = reshape([s.name for s in spectra], (1, length(spectra)))
    plot(shifts, [resample(s, shifts) for s in spectra]; lab=labels)
    plot!(xflip=true, yticks=[], grid=false, yforeground_color_axis=false; kw...)
end

function integral_curve(a::AbstractArray{T,1}, scale::Float64; shift=0.0) where T
    x = similar(a)
    x[1] = shift
    for i in eachindex(a[2:end])
        x[i+1] = x[i] + a[i]*scale
    end
    return x
end

function integration_plot(s::Spectrum, labels::AbstractArray{T}) where T <: AbstractString
    ints = integrate(s)
    if length(ints)==0
        return
    end
    maxint = maximum(ints)
    maxpoint = maximum(maximum.(intrng_data(s)))
    scale = maxpoint/maxint
    shifts = intrng_shifts(s)
    plot!(shifts, [integral_curve(s[ppmtoindex(s,r)], scale) for r in s[s.default_proc].intrng], color=:red, linewidth=1.5, legend=false)
    annotation_height = -0.2maxpoint
    for (i,r) in enumerate(s[s.default_proc].intrng)
        annotate!([(mean(r), annotation_height, text(labels[i], font(10,90.0,colorant"red")))])
    end
end

function integration_plot(s::Spectrum)
    n = length(s[s.default_proc].intrng)
    labels = collect(string(c) for c in 'a':('a'+n-1))
    integration_plot(s, labels)
end