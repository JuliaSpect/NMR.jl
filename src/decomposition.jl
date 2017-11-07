"""
    overlay!(signal, chunks, positions)

Overlay a series of chunks on top of a signal at the
specified positions.
"""
function overlay!(signal, chunks, positions)
    for (i,c) in enumerate(chunks)
        p = positions[i]
        signal[p:(p+length(c)-1)] .+= c
    end
    signal
end


function candidates(signal, chunk, start_pos; tol = 125)
    l = length(chunk)
    positions = (-tol+start_pos):(tol+start_pos)
    s = normalize(signal[(-tol+start_pos):(tol+start_pos+l-1)])
    c = normalize(chunk)
    corr = zeros(tol*2+1)
    for i in eachindex(corr)
        corr[i] = sum(c.*s[i:(i+l-1)])
    end
    cl = length(corr)
    [positions[i] for i in eachindex(positions) 
        if corr[i]>0.2 && (i==cl || corr[i] > corr[i+1]) && 
                          (i==1 || corr[i] > corr[i-1])]
end
 
candidates(s::Spectrum, args...; kw...) = candidates(s[:],args...;kw...)
candidates(s::Spectrum, l::Spectrum, args...; kw...) =
    [ candidates(s[:], l[r], r.start, args...; kw...) for r in intrng_indices(l) ]


function cand_signals(s::Spectrum, l::Spectrum)
    positions = candidates(s, l)
    chunks = [l[r] for r in intrng_indices(l)]
    (overlay!(zeros(length(s[:])), chunks, comb) for
        comb in Base.product(positions...))
end

function guess_matrices(s::Spectrum, lib::Array{Spectrum,1})
    gens = [cand_signals(s,l) for l in lib]
    indices = [i for i in eachindex(gens) if length(gens[i])>0]
    product_gens = Base.product([g for g in gens if length(g)>0]...)
    (indices, (hcat(p...) for p in product_gens))
end

function lsq_analyze(s::Vector{Float64}, l::Matrix{Float64})
    soln = l\s # Least-squares solution
    recon = l*soln
    res = norm(recon .- s)/norm(s)
    (soln,res)
end

struct DecompositionResult
    coefficients :: Vector{Float64}
    refnums :: Vector{Int}
    signal :: Vector{Float64}
    matrix :: Matrix{Float64}
end

function lsq_analyze(s::Spectrum, lib::Array{Spectrum,1}, dark_areas::Vector{Tuple{Float64,Float64}} = Tuple{Float64,Float64}[])
    sig = copy(s[:])
    for a in dark_areas
        r = ppmtoindex(s,a[1]):ppmtoindex(s,a[2])
        sig[r] = 0.0
    end
    refnums, matrices = guess_matrices(s, lib)
    res = vec(collect(lsq_analyze(sig, m) for m in matrices))
    _,best = findmin(r[2] for r in res)
    for (i,m) in enumerate(matrices)
        if i==best
            return DecompositionResult(res[best][1], refnums, sig, m)
        end
    end
end
