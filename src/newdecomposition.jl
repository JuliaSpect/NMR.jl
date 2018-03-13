const MINFACT = 0.05
const HITCORR = 0.9
const MINCORR = 0.3 # reasonable default
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

function alignments(signal, chunk, start_pos; tol=250, fuzziness=1.5,
                    minsig=MINFACT*mean(signal), mincorr=MINCORR)
    l = length(chunk)
    section = max(1, -tol+start_pos):min(length(signal), tol+l+start_pos)
    sig = @view signal[section]
    # reject weak signals, important if signal only contains noise
    if mean(sig) < minsig
        return Tuple{Int64,Float64}[]
    end
    s = normalize(sig)
    c = normalize(chunk)
    corr = zeros(length(section)-l+1)
    for i in eachindex(corr)
        corr[i] = sum(c.*s[i:(i+l-1)])
    end
    M = maximum(corr)
    ps = section.start - 1
    Tuple{Int64,Float64}[ (i+ps, corr[i]) for i=2:length(corr) if corr[i] > mincorr &&
                                           corr[i] > corr[i-1] &&
                                           corr[i] > corr[i+1] &&
                                           M/corr[i] > fuzziness ]
end

alignments(s::Spectrum, args...; kw...) = alignments(s[:],args...;kw...)
alignments(s::Spectrum, l::Spectrum, args...; kw...) =
    [ alignments(s[:], l[r], r.start, args...; kw...) for r in intrng_indices(l) ]


function aligned_signals(s::Spectrum, l::Spectrum; sloppiness=0, kw...)
    positions = alignments(s, l; kw...)
    matched = filter(p->!isempty(p), positions)
    if length(matched) == 0 ||  length(positions) - length(matched) > sloppiness
        return []
    end
    chunks = [l[r] for (i,r) in enumerate(intrng_indices(l)) if !isempty(positions[i])]
    (overlay!(zeros(length(s[:])), chunks, comb) for
        comb in Base.product(matched...))
end

function guess_matrices(s::Spectrum, lib::Array{Spectrum,1}; exclude = [], minfact=fill(MINFACT, length(lib)), sloppiness=Dict(), kw...)
    # println("""Ratio factors: $(join(minfact, ", "))""")
    sk = keys(sloppiness)
    gens = [cand_signals(s,l; minfact=minfact[i], sloppiness=get(sloppiness,i,0), kw...) for (i,l) in enumerate(lib)]
    indices = [i for i in eachindex(gens) if length(gens[i])>0 && i ∉ exclude]
    if isempty(indices)
        # no match found
        return ([], [])
    end
    # product_gens = Base.product([g for g in gens if length(g)>0]...)
    product_gens = Base.product(gens[indices]...)
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

function analyze(s::Spectrum, lib::Array{Spectrum,1}; kw...)
    found = []
    signal = s[:]

    lsq_analyze(s, lib, found)
end

function lsq_analyze(s::Spectrum, lib::Array{Spectrum,1}; kw...)
    found = Int[]
    coeffs = Float64[]
    mincorrs = CORR_THRESHOLDS[:]
    mincorr = shift!(mincorrs)
    vecs = Float64[]
    ss = deepcopy(s)
    sig = ss[:]
    ms = mean(sig)
    println("Signal mean: $ms")
    while true
        refnums, matrices = guess_matrices(ss, lib; avg=ms, mincorr = mincorr, exclude=found, kw...)
        println("Forming $(length(matrices)) matrices (correlation threshold = $mincorr)")
        if issubset(refnums, found)
            if isempty(mincorrs)
                return DecompositionResult(coeffs,found,s[:],
                                           isempty(vecs)?Matrix{Float64}(0,0):reshape(vecs,(:,length(found))))
            else
                mincorr = shift!(mincorrs)
                continue
            end
        end
        res = vec(collect(lsq_analyze(sig, m) for m in matrices))
        _,best = findmin(r[2] for r in res)
        # _,max_component = findmax(res[best][1]) # Largest constituent
        for (i,m) in enumerate(matrices)
            if i==best
                # return DecompositionResult(res[best][1], refnums, sig, m)
                for (j,r) in enumerate(refnums)
                # r = refnums[max_component]
                    # println("j: $j, max_component: $max_component")
                    sig .-= m[:,j] .* res[i][1][j]
                    println("Found: #$r")
                    if r ∉ found
                        append!(coeffs, res[i][1][j])
                        append!(found, r)
                        append!(vecs, m[:,j])
                    end
                end
                break
                # end
            end
        end
        ss[:] = sig
    end
end

function decompose(d::DecompositionResult)
    if isempty(d.coefficients)
        return ([],zeros(length(d.signal)), d.signal)
    end
    recon = d.matrix * d.coefficients
    residue = d.signal .- recon
    components = d.coefficients'.*d.matrix
    (components, recon, residue)
end
