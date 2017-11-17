module D

using NMR
using NMRquant
using IterTools
import NMR: intrng_indices

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
function candidates(s::Spectrum, lib::Array{Spectrum,1}, args...; kw...)
    cands = [ collect(product(candidates(s, l, args...; kw...)...)) for l in lib ]
    collect(product(c for c in cands if length(c)>0))
end

# function cand_signals(s::Spectrum, l::Spectrum)
#     positions = candidates(s, l)
#     chunks = [l[r] for r in intrng_indices(l)]
#     (overlay!(zeros(length(s[:])), chunks, comb) for
#         comb in product(positions...))
# end

function guess_combinations(s::Spectrum, lib::Array{Spectrum,1})
    gens = [product(candidates(s,l)...) for l in lib]
    indices = [i for i in eachindex(gens) if length(gens[i])>0]
    product_gens = product([g for g in gens if length(g)>0]...)
    # result = zeros(length(s), length(gens))
    (indices, product_gens)
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
    l = length(s)
    for a in dark_areas
        r = NMR.ppmtoindex(s,a[1]):NMR.ppmtoindex(s,a[2])
        sig[r] = 0.0
    end
    refnums, matrices = guess_combinations(s, lib)
    # res = collect(lsq_analyze(sig, m) for m in matrices)
    res = pmap(guess -> begin
        x = zeros(l,length(refnums))
        for (n,g) in zip(refnums,guess)
            overlay!(x, g, NMR.intrng_data(lib[n]))
        end
        lsq_analyze(sig, )
    end)
    _,best = findmin(r[2] for r in res)
    for (i,m) in enumerate(matrices)
        if i==best
            return DecompositionResult(res[best][1], refnums, sig, m)
        end
    end
end

end