"""
    baseline_correct!(s::Spectrum, ppm_range::Tuple{Float64,Float64})

Subtract the mean of a defined region `ppm_range` from the Spectrum `s`
without altering the overall shape of the baseline.
"""
function baseline_correct!(s::Spectrum, ppm_range::Tuple{Float64,Float64})
    s[:] .-= mean(s[ppm_range])
    s
end
