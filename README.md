# NMR
This package helps with analysis and visualization of NMR data.

![screenshot](screenshot.png)

## Getting started

Make sure you have Julia 0.6 or better installed.

```julia
Pkg.init()  # if you are installing a package for the first time
Pkg.clone("https://github.com/Julia-NMR/NMR.jl")

using NMR
s = Spectrum("<path to Bruker format experiment>/1")
s.acqu   # Acquisition parameters
s[(2.3,1.5)]   # Data points in the 1.5–2.3 ppm range
plot(s)
```

## Implemented features
* Reading and writing Bruker FIDs and processed data
* Importing of acquisition and processing parameters and integration ranges
* Signal interpolation
* Arithmetic operations between NMR spectra
* Signal integration and integral display
* Utility functions for plotting with
  [Plots.jl](https://github.com/JuliaPlots/Plots.jl)
* Least squares signal decomposition

## TODO
* Peak detection
* Window functions
* Automatic phase correction

## Authors
* [S. Hessam M. Mehr](https://hessammehr.github.io) ([github](https://github.com/hessammehr))

## License
MIT

Contributions are very welcome
