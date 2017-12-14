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


function candidates(signal, chunk, start_pos; tol = 200)
    l = length(chunk)
    positions = (-tol+start_pos):(tol+start_pos)
    s = vcat(zeros(tol),signal[start_pos:start_pos+l-1],zeros(tol))
    # reject weak signals, important if signal only contains noise
    if norm(s)/norm(chunk) < 0.05
        return []
    end
    normalize!(s)
    c = normalize(chunk)
    corr = zeros(tol*2+1)
    for i in eachindex(corr)
        corr[i] = sum(c.*s[i:(i+l-1)])
    end
    cl = length(corr)
    # good debug hook
    # Gaussian window to favor smaller shifts
    corr .*= exp.(-2(-tol:tol).^2./tol^2)
    [positions[i] for i in 2:(cl-1)
        if corr[i]>0.3 && (corr[i] > corr[i+1]) && (corr[i] > corr[i-1])]
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
    if sum(length(g) for g in gens) == 0
        # no match found
        return ([], [])
    end
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

function lsq_analyze(s::Spectrum, lib::Array{Spectrum,1})
    found = Int[]
    coeffs = Float64[]
    vecs = Float64[]
    ss = deepcopy(s)
    sig = ss[:]
    # sig = s[:]
    # for a in dark_areas
    #     r = ppmtoindex(s,a[1]):ppmtoindex(s,a[2])
    #     sig[r] = 0.0
    # end
    while true
        refnums, matrices = guess_matrices(ss, lib)
        if issubset(refnums, found)
            return DecompositionResult(coeffs,found,s[:],
                                       isempty(vecs)?Matrix{Float64}(0,0):reshape(vecs,(:,length(found))))
        end
        res = vec(collect(lsq_analyze(sig, m) for m in matrices))
        _,best = findmin(r[2] for r in res)
        for (i,m) in enumerate(matrices)
            if i==best
                # return DecompositionResult(res[best][1], refnums, sig, m)
                for (j,r) in enumerate(refnums)
                    if r ∉ found
                        sig .-= m[:,j] .* res[i][1][j]
                        append!(coeffs, res[i][1][j])
                        append!(found, r)
                        append!(vecs, m[:,j])
                        break
                    end
                end
            end
        end
        ss[:] = sig
    end
end

function decompose(d::DecompositionResult)
    recon = d.matrix * d.coefficients
    residue = d.signal .- recon
    components = d.coefficients'.*d.matrix
    (components, recon, residue)
end
