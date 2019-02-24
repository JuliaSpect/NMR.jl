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

@recipe function f(s::Union{Spectrum, AbstractArray{Spectrum}};
                   integrate=false, Δ=union_range(s), npoints=0)
    isarray = isa(s, AbstractArray)
    ss = isarray ? s : [s]
    npoints = npoints == 0 ? length(ss[1]) : npoints
    xflip --> true
    grid --> false
    legend --> isarray
    framestyle --> :box
    xticks --> :auto
    yticks --> []
    xlabel --> "Chemical shift (ppm)"
    for s in ss
        @series begin
            x = range(Δ[1]; stop=Δ[2], length=npoints)
            y = resample(s, x)
            x,y
        end
        if integrate
            @series begin
                primary := false
                color --> :red
                legend := false
                linewidth --> 1.5
                integral_curve(s,Δ)
            end
        end
    end
end

function integral_curve(a::AbstractArray{T,1}, scale::Float64; shift=0.0) where T
    isempty(a) && return Float64[]
    x = similar(a)
    x[1] = shift
    for i in eachindex(a[2:end])
        x[i+1] = x[i] + a[i]*scale
    end
    return x
end

function integral_curve(s::Spectrum, Δ=limits(s))
    ints = integrate(s)
    if length(ints)==0
        return
    end
    maxint = maximum(ints)
    maxpoint = maximum(maximum.(intrng_data(s)))
    scale = maxpoint/maxint
    shifts = filter.(δ->minimum(Δ) <= δ <= maximum(Δ), intrng_shifts(s))
    (shifts, [integral_curve(s[ppmtoindex(s,r)], scale) for r in shifts])
end