import Plots: plot, plot!, @colorant_str, annotate!, text, font
using Formatting

plot(s::Spectrum; integrate = false, kw...) = begin
    lo,hi = limits(s)
    shifts = linspace(hi, lo, length(s))
    p = plot(color=colorant"#0a3faa",shifts, s[:])
    if integrate == true
        integration_plot(s, 1)
    end
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

function integration_plot(s::Spectrum, refindex::Int)
    ints = integrate(s)
    maxint = maximum(ints)
    maxpoint = maximum(maximum.(intrng_data(s)))
    scale = maxpoint/maxint
    shifts = intrng_shifts(s)
    plot!(shifts, [integral_curve(s[ppmtoindex(s,r)], scale) for r in s[s.default_proc].intrng], color=:red, linewidth=1.5, legend=false)
    annotation_height = -0.2maxpoint
    for (i,r) in enumerate(s[s.default_proc].intrng)
        annotate!([(mean(r),annotation_height,text(sprintf1("%.2f",ints[i]/ints[refindex]),font(10,90.0,colorant"red")))])
    end
end