import StatsBase

const MINFACT = 0.1
const HITCORR = 0.9
const MINCORR = 0.25 # reasonable default
const MAXFUZZ = 1.8
const MINFUZZ = 1.0
const TARGETNGUESS = 200
const Guess{N} = NTuple{N,Tuple{Int64,Float64}} where N

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
                    minsig=MINFACT*norm(signal), mincorr=MINCORR)
    l = length(chunk)
    section = max(1, -tol+start_pos):min(length(signal), tol+l-1+start_pos)
    s = @view signal[section]
    # reject weak signals, important if signal only contains noise
    c = normalize(chunk)
    # return s,c
    corr = zeros(length(section)-l+1)
    lc = length(corr)
    lratio = l/length(signal)
    for i in eachindex(corr)
        sig = @view s[i:(i+l-1)]
        xc = dot(c,sig)
        corr[i] = xc < minsig*lratio ? 0.001 :  xc/norm(sig)
    end
    M = maximum(corr)
    if M < mincorr
        return Tuple{Int64,Float64}[]
    end
    # println(M)
    ps = section.start - 1
    let lc=lc, corr=corr, ps=ps, M=M
        Tuple{Int64,Float64}[ (i+ps, corr[i]) for i=eachindex(corr) if
                              # corr[i] > mincorr &&
                              (i==1 || corr[i] > corr[i-1]) &&
                              (i==lc || corr[i] > corr[i+1]) &&
                              M/corr[i] < fuzziness ]
    end
end

alignments(s::Spectrum, args...; kw...) = alignments(s[:],args...;kw...)
alignments(s::Spectrum, l::Spectrum, args...; kw...) =
    Array{Tuple{Int64,Float64},1}[ alignments(s[:], l[r], r.start, args...; kw...) for r in intrng_indices(l) ]

guesses(s::Spectrum, l::Spectrum; kw...) = vec(collect(Base.product(alignments(s, l; kw...)...)))
function guesses_adaptive(s::Spectrum, l::Spectrum;
                          ntarget=TARGETNGUESS, δ=ntarget/4, kw...)
    fuzz = (MAXFUZZ + MINFUZZ) / 2
    @binary_opt ( *(length.(alignments(s, l; fuzziness=fuzz))...) - ntarget ) fuzz MINFUZZ MAXFUZZ δ
    guesses(s, l; fuzziness=fuzz)
end

function projection(signal, chunk, start_pos)
    l = length(chunk)
    s = @view signal[start_pos:start_pos+l-1]
    dot(s, chunk)/norm(chunk)^2
end

positions(guess) = Int64[ g[1] for g in guess ]
fitness(guess) = Float64[ g[2] for g in guess ]

function fitness_weights(guess, λ=10.0)
    f = fitness(guess)
    F = maximum(f)
    exp.(λ/F*(f.-F))
end

function projection_weights(projs, fitness_weights=ones(length(projs)), η=2.0)
    # p = projections(s, l, positions(guess))
    fw = StatsBase.weights(fitness_weights)
    σ = std(projs, fw; corrected=false)/η
    m = mean(projs, fw)
    exp.(-(projs.-m).^2/2σ^2)
end

function projections(s::Spectrum, l::Spectrum, positions)
    sig = s[:]
    Float64[ projection(sig, r, p) for (p,r) in zip(positions, intrng_data(l)) ]
end

# function projection(s::Spectrum, l::Spectrum, positions, weights=ones(positions))
#     projs = projections(s, l, positions)
#     mean(projs, StatsBase.weights(weights))
# end

function projection(s::Spectrum, l::Spectrum, guess::NTuple{N,Tuple{Int64,Float64}}) where N
    projs = projections(s, l, positions(guess))
    if length(projs) == 1
        return first(projs)
    end
    fw = fitness_weights(guess)
    # projection(s, l, positions(guess), f(s, l, guess))
    mean(projs, StatsBase.weights(projection_weights(projs, fw)))
    # projection(s, l, positions(guess), Float64[f(ff) for ff in fitness(guess)])
end

function projection_score(s::Spectrum, l::Spectrum, guess::NTuple{N,Tuple{Int64,Float64}}) where N
    projs = projections(s, l, positions(guess))
    mean(projs)/std(projs, StatsBase.weights(fitness(guess)); corrected = false)
end

function fit_score(s::Spectrum, l::Spectrum, guess::NTuple{N,Tuple{Int64,Float64}}) where N
    if N == 0
        0.0
    else
        *((g[2] for g in guess)...)^(1/N)
    end
end

# overall score, not comparable between references
score(s, l, g) = fit_score(s, l, g) * projection_score(s, l, g)

synthesize(l::Spectrum, positions) =
    overlay!(zeros(length(l)), intrng_data(l), positions)
synthesize(l::Spectrum, guess::Guess{N}) where N =
    synthesize(l, positions(guess))

# function score(s::Spectrum, l::Spectrum, guess::AbstractArray{Tuple{Int,Float64}})
#     *((g[2] for g in guess)...)
# end

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

struct DecompositionResult
    coefficients :: Vector{Float64}
    refnums :: Vector{Int}
    signal :: Vector{Float64}
    matrix :: Matrix{Float64}
end

function lsq_analyze(s::Spectrum, lib::AbstractArray{Spectrum}, found; kw...)
    gs = guesses_adaptive.(s, lib; kw...)
    fit_scores = let s=s
        Array{Float64,1}[ isempty(g) || i ∈ found ? [0.0] : fit_score.(s, l, g)
                          for (i,l,g) in zip(1:length(lib),lib, gs) ]
    end
    scores = let s=s
        Array{Float64,1}[ isempty(g) || i ∈ found ? [0.0] : score.(s, l, g)
                          for (i,l,g) in zip(1:length(lib),lib, gs) ]
    end
    # find best guess per reference based on overall score
    bestinds = indmax.(scores)
    # find best reference based on fit_score
    max_score,max_ref = findmax(getindex.(fit_scores, bestinds))
    if max_score == 0.0
        return (s, 0, (), 0.0, Float64[])
    end
    # pss = score.(s, lib[maxind], gs[maxind])
    # _,maxguessind = findmax(pss)
    maxguess = gs[max_ref][bestinds[max_ref]]
    ss = deepcopy(s)
    l = synthesize(lib[max_ref], maxguess)
    p = projection(s, lib[max_ref], maxguess)
    ss[:] .-= p.*l
    ss, max_ref, maxguess, p, l.*p
end

function lsq_analyze(s::Spectrum, lib::AbstractArray{Spectrum};
                     callback = refs -> println("Found: #$refs"), kw...)
    found = Int64[]
    coeffs = Float64[]
    vecs = Array{Float64,1}[]
    ss = copy(s[:])
    while true
        s, m, g, p, v = lsq_analyze(s, lib, found; kw...)
        if m == 0
            break
        end
        push!(found, m)
        push!(vecs, v)
        push!(coeffs, p)
        callback([m])
    end
    DecompositionResult(coeffs, found, ss, hcat(vecs...))
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
