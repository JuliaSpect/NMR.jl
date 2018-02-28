import Plots: plot, plot!, @colorant_str, annotate!, text, font
using Formatting
using RecipesBase

function plot_limits(s::NMR.Spectrum, margin = 0.1)
    m, M = (minimum(minimum(d) for d in intrng_data(s)),
            maximum(maximum(d) for d in intrng_data(s)))
    Δ = M - m
    (m - margin*Δ, M + margin * Δ)
end

function plot_limits(a::AbstractArray, margin = 0.1)
    m, M = extrema(a)
    Δ = M - m
    (m - margin*Δ, M + margin * Δ)
end

# plot(s::Spectrum; integrate = false, kw...) = begin
#     lo,hi = limits(s)
#     shifts = linspace(hi, lo, length(s))
#     p = plot(shifts, s[:], framestyle=:box)
#     if integrate == true
#         integration_plot(s)
#     end
#     plot!(xflip=true, yticks=[], leg=false; kw...)
# end

@recipe function f(s::Spectrum; integrate=false, δ=limits(s), npoints = length(s))
    xflip --> true
    grid --> false
    legend --> false
    framestyle --> :box
    yticks --> []
    xlims := δ
    @series begin
        x = linspace(δ[1], δ[2], npoints)
        y = resample(s, x)
        x,y
    end
    if integrate
        @series begin
            primary := false
            color --> :red
            legend := false
            linewidth --> 1.5
            integral_curve(s)
        end
    end
end

# function plot(spectra::AbstractArray{Spectrum}; kw...)
#     shifts = union_shifts(spectra)
#     labels = reshape([s.name for s in spectra], (1, length(spectra)))
#     plot(shifts, [resample(s, shifts) for s in spectra]; lab=labels)
#     plot!(xflip=true, yticks=[], grid=false, yforeground_color_axis=false; kw...)
# end

function integral_curve(a::AbstractArray{T,1}, scale::Float64; shift=0.0) where T
    x = similar(a)
    x[1] = shift
    for i in eachindex(a[2:end])
        x[i+1] = x[i] + a[i]*scale
    end
    return x
end

function integral_curve(s::Spectrum)
    ints = integrate(s)
    if length(ints)==0
        return
    end
    maxint = maximum(ints)
    maxpoint = maximum(maximum.(intrng_data(s)))
    scale = maxpoint/maxint
    shifts = intrng_shifts(s)
    (shifts, [integral_curve(s[ppmtoindex(s,r)], scale) for r in intrng(s)])
    # annotation_height = -0.2maxpoint
    # for (i,r) in enumerate(intrng(s))
    #     annotate!([(mean(r), annotation_height, text(labels[i], font(10,90.0,colorant"red")))])
    # end
end

# function integration_plot(s::Spectrum)
#     n = length(intrng(s))
#     labels = collect(string(c) for c in 'a':('a'+n-1))
#     integration_plot(s, labels)
# end