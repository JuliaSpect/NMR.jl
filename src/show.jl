import Base: show

function show(io::IO, s::Spectrum)
    println(io, "==== $(s.name) ($(s.expno)) ====")
    println(io, "FID size: $(length(s.fid))")
    println(io, "# processings: $(length(s.procs))") 
    println(io, "Default proc. no.: $(s.default_proc)")
    out = PipeBuffer()
    for k in sort(collect(keys(s.procs)))
        show(out, s.procs[k])
        println(out, "")
    end
    for l in eachline(out)
        println(io, "    $l")
    end
end

function show(io::IO, p::ProcessedSpectrum)
    println(io, "---- Proc. no. $(p.procno) ----")
    println(io, strip(p.title))
    println(io, "Processed size: $(p["SI"]) points")
end