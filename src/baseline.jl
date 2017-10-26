function baseline_correct!(s::Spectrum, ppm_range::Tuple{Float64,Float64})
    s[:] .-= mean(s[ppm_range])
    s
end
