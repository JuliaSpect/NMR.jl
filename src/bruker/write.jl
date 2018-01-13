import Base: show

function show(io::IO, intrng::Array{Intrng,1})
    write(io, """A 1.0 #regions in PPM
# low field   high field  bias        slope\n""")
    entries = sort(intrng, by=(Δ->Δ[1]), rev=true)
    for (i,(δ1, δ2)) in enumerate(entries)
        write(io, " $δ1  $δ2  0.0  0.0  # for region $i\n")
    end
end

function dump(path, s::Spectrum)
    if !isdir(path)
        mkdir(path)
    end

    # Save each child proc
    for (n,p) in s.procs
        dump(joinpath(path, "pdata", string(n)), p)
    end

    # Save fid
    write(open(joinpath(path,"fid"), "w"), Int32.(round.(s.fid)))

    # TODO: Save acqu/acqus
end

function dump(path, s::Spectrum, templatepath)
    cp(templatepath, path)
    dump(path, s)
end

function dump(path, p::ProcessedSpectrum)
    if !isdir(path)
        mkdir(path)
    end

    # Save 1r and 1i
    open(joinpath(path, "1r"), "w") do f
        write(f, Int32.(round.(p.re_ft)))
    end
    open(joinpath(path, "1i"), "w") do f
        write(f, Int32.(round.(p.im_ft)))
    end

    # Save intrng
    open(joinpath(path, "intrng"), "w") do f
        show(f, p.intrng)
    end

    # Save title
    open(joinpath(path, "title"), "w") do f
        write(f, p.title)
    end

    # TODO: Save proc/procs
end