# Gibbs sampler
# -------------
abstract type MHProposal end

function lognormalize(l)
    x = exp.(l .- maximum(l))
    return x ./ sum(x)
end

struct BetaProposal{T}
    d::Beta{T}
    BetaProposal(α::T) where T = new{T}(Beta(α,α))
end

function propose(d::BetaProposal, x)
    xp = rand(d.d)
    qp = logpdf(d.d, x) - logpdf(d.d, xp)
    xp, qp
end

struct BetaSwitchProposal{T}
    d::Beta{T}
    r::Float64  # switch probability
    BetaSwitchProposal(α::T, r) where T = new{T}(Beta(α,α), r)
end

function propose(d::BetaSwitchProposal, x::Real)
    xp = rand(d.d)
    qp = logpdf(d.d, x) - logpdf(d.d, xp)
    xp, qp
end

struct BetaFlipProposal{T}
    d::Beta{T}
    r::Float64  # switch probability
    rj::Float64
    BetaFlipProposal(α::T, r, rj) where T = new{T}(Beta(α,α), r, rj)
end

function propose(d::BetaFlipProposal, x::Real)
    if rand() < d.rj
        return 1-x, 0.
    else
        xp = rand(d.d)
        qp = logpdf(d.d, x) - logpdf(d.d, xp)
        return xp, qp
    end
end

haldane(y) = 0.5*(1-exp(-y))
invhaldane(x) = -log(1 - 2x)

function rrates(xs, j)
    ys = invhaldane.(xs)
    left = reverse(cumsum(ys[j-1:-1:1]))
    rght = cumsum(ys[j:end-1])
    [haldane.(left); NaN; haldane.(rght)]
end

function gff(A::Architecture, p, y, j)
    @unpack loci, rrate = A
    rs = rrates(rrate, j)
    x = 0.
    for i=1:length(p)
        i == j && continue
        @unpack s1, s01, s11 = loci[i]
        sa = s1 + s01
        sb = s11 - 2s01
        x += (y[i]-(1-p[i]))*(sa + sb*(1-p[i])) / rs[i]
    end
    return exp(x)  # XXX factor 2 disappears here with linkage -> in the rᵢ
end

function unnormϕ(M, ps, j)
    g = gff(M.arch, ps, M.y, j)  # gene flow factor due to L-1 other loci
    return unnormϕ(M, M.arch.loci[j], ps[j], M.y[j], g)
end

function unnormϕ(M::MainlandIslandModel, locus::HapDipLocus, p, y, g) 
    @unpack m, u = M
    @unpack s1, s01, s11 = locus
    Ne = getNe(M)
    A  = 2Ne*(u + m*g*y) - 1
    B  = 2Ne*(u + m*g*(1-y)) - 1
    sa = s1 + s01
    sb = s11 - 2s01
    fp = Ne*p*(2sa + sb*(2 - p))
    return B*log(p) + A*log(1-p) - fp 
end

function gibbs(model, kernel, p::Vector{T}, n; drop=0) where T
    L = length(p)
    ℓ = T[unnormϕ(model, p, j) for j=1:L]
    P = Matrix{T}(undef, n, L)
    Λ = Vector{T}(undef, n)
    for i=1:n 
        p, ℓ = joint_update(model, kernel, p, ℓ) 
        for j=1:length(p)
            p, ℓ = gibbs_update(model, kernel, p, ℓ, j)
        end
        P[i,:] = p
        Λ[i]   = sum(ℓ) 
    end
    1 .- P[drop+1:end,:], Λ  
    # We are sampling the frequency of the 0 allele (typically locally
    # beneficial). For consistency with individual-based simulations we should
    # report the frequency of the 1 allele (usually referred to as q)
end

function gibbs_update(model, kernel, p, ℓ, j)
    _p = copy(p)
    _ℓ = copy(ℓ)
    pj, qj = propose(kernel, p[j])
    _p[j] = pj
    _ℓ[j] = unnormϕ(model, _p, j)
    if !isfinite(_ℓ[j]) || !(log(rand()) < _ℓ[j] - ℓ[j] + qj)  # reject
        _p[j] = p[j]
        _ℓ[j] = ℓ[j]
    end
    return _p, _ℓ
end

joint_update(model, kernel::BetaProposal, p, ℓ) = p, ℓ

function joint_update(model, kernel, p, ℓ)
    (rand() > kernel.r) && return (p, ℓ) 
    _p = 1 .- p
    _ℓ = [unnormϕ(model, _p, j) for j=1:length(p)]
    if !(log(rand()) < sum(_ℓ) - sum(ℓ))  # reject
        _p = p
        _ℓ = ℓ
    end
    return _p, _ℓ
end

# -----------------------------
# New take using adaptive stuff

struct GibbsSampler{P}
    proposals::Vector{P}
end

Base.show(io::IO, g::GibbsSampler) = write(io, "GibbsSampler(L=$(length(g.proposals)))")

function gibbs(model, smpler::GibbsSampler, p::Vector{T}, n; drop=0) where T
    L = length(p)
    P = Matrix{T}(undef, n, L)
    for i=1:n 
        p = joint_update(model, smpler, p)
        for j=1:length(p)
            p = gibbs_update(model, smpler, p, j)
            p = switch_update(model, smpler, p, j)
        end
        P[i,:] = p
    end
    1 .- P[drop+1:end,:]
    # We are sampling the frequency of the 0 allele (typically locally
    # beneficial). For consistency with individual-based simulations we should
    # report the frequency of the 1 allele (usually referred to as q)
end

function gibbs_update(model, smpler::GibbsSampler, p, j)
    # compute the present factor (before proposal -- XXX the factor has changed
    # since we last saw it due to the updates at other sites!)
    ℓ = unnormϕ(model, p, j)
    # make a copy of the p vector and do a proposal at site j
    prop = smpler.proposals[j]
    _p = copy(p)
    pj, qj = prop(p[j])
    _p[j] = pj
    _ℓ = unnormϕ(model, _p, j)
    if !isfinite(_ℓ) || !(log(rand()) < _ℓ - ℓ + qj)  # reject
        _p[j] = p[j]
    else
        accept!(prop)
    end
    return _p
end

function switch_update(model, smpler::GibbsSampler, p, j)
    ℓ = unnormϕ(model, p, j)
    _p = copy(p)
    _p[j] = 1 - p[j]
    _ℓ = unnormϕ(model, _p, j)
    if !isfinite(_ℓ) || !(log(rand()) < _ℓ - ℓ)  # reject
        _p[j] = p[j]
    end
    return _p
end

function joint_update(model, smpler::GibbsSampler, p)
   _p = 1 .- p
    ℓ = [unnormϕ(model,  p, j) for j=1:length(p)]
   _ℓ = [unnormϕ(model, _p, j) for j=1:length(p)]
    L = sum( ℓ)
   _L = sum(_ℓ)
    if !isfinite(_L) || !(log(rand()) < _L - L)  # reject
        _p = p
    end
    return _p
end

