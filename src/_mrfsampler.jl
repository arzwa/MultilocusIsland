# Gibbs sampler for the finite island model

# Finite island model MRF approach, unnormalized density at locus j in deme k
function ϕ(M::FiniteIslandModel, p, j, k)
    d = ndeme(M)
    y = zeros(nloci(M))
    for i=1:d
        y .+= p[i,:] * M.m[i]  # contribution to the gamete pool
    end
    y ./= sum(M.m)
    # normalize, this will be the allele frequency in the migrant pool (what
    # about the heterozygosity)? Migration at diploid level? ...
    m = M.m[k]
    @unpack Ne, u, A = M.D[k]
    @unpack s1, s01, s11 = A[j]
    g  = gff(A, p[k,:], p[k,:] .* (1 .- p[k,:]), y, j)
    A  = 2Ne*(u + m*g*y[j]) - 1
    B  = 2Ne*(u + m*g*(1-y[j])) - 1
    sa = s1 + s01
    sb = s11 - 2s01
    fp = Ne*p[k,j]*(2sa + sb*(2 - p[k,j]))
    return B*log(p[k,j]) + A*log(1-p[k,j]) - fp 
end

struct GibbsSampler{P}
    proposals::Matrix{P}
end

function GibbsSampler(M::FiniteIslandModel)
    GibbsSampler([UnitIntervalProposal() for k=1:ndeme(M), i=1:nloci(M)])
end

Base.show(io::IO, g::GibbsSampler) = write(io, "GibbsSampler(L=$(length(g.proposals)))")

function gibbs(model, smpler::GibbsSampler, p::Matrix{T}, n; drop=0) where T
    d, L = size(p)
    P = Array{T,3}(undef, n, d, L)
    @showprogress 0.1 "[Sampling]" for i=1:n 
        for j=1:L, k=1:d
            p = gibbs_update(model, smpler, p, j, k)
        end
        P[i,:,:] = p
    end
    1 .- P[drop+1:end,:,:]
    # We are sampling the frequency of the 0 allele (typically locally
    # beneficial). For consistency with individual-based simulations we should
    # report the frequency of the 1 allele (usually referred to as q)
end

function gibbs_update(model, smpler::GibbsSampler, p, j, k)
    # compute the present factor (before proposal -- XXX the factor has changed
    # since we last saw it due to the updates at other sites!)
    ℓ = ϕ(model, p, j, k)
    # make a copy of the p vector and do a proposal at site j
    prop = smpler.proposals[k,j]
    _p = copy(p)
    pj, qj = prop(p[k,j])
    _p[k,j] = pj
    _ℓ = ϕ(model, _p, j, k)
    if !isfinite(_ℓ) || !(log(rand()) < _ℓ - ℓ + qj)  # reject
        _p[k,j] = p[k,j]
    else
        accept!(prop)
    end
    return _p
end

