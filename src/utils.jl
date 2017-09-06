limits(s :: Spectrum) = limits(s.acqu["O1P"], s.acqu["SW"])
limits(o1p, sw) = (o1p - sw/2, o1p + sw/2)

### Plotting functions

import Plots: plot

plot(s :: Spectrum, fid :: Bool = false) = begin
    if fid
    lo,hi = limits(s)
    shifts = linspace(lo, hi, length(s.re_ft))    
    # (shifts, s.re_ft)
    plot(shifts, s.re_ft, xflip=true)
end
