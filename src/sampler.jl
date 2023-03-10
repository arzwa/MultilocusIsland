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

function gff(loci::Vector{<:HapDipLocus}, p, y, j)
    x = 0.
    for i=1:length(p)
        i == j && continue
        @unpack s1, s01, s11 = loci[i]
        x += (y[i]-(1-p[i]))*(s1 + s01 + (s11 - 2s01)*(1-p[i]))
    end
    return exp(2x)
end

function unnormϕ(M, ps, j)
    g = gff(M.arch, ps, M.y, j)  # gene flow factor due to L-1 other loci
    return unnormϕ(M, M.arch[j], ps[j], M.y[j], g)
end

function unnormϕ(M::MainlandIslandModel, locus::HapDipLocus, p, y, g) 
    @unpack Ne, m, u = M
    @unpack s1, s01, s11 = locus
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
    if !(log(rand()) < _ℓ[j] - ℓ[j] + qj)  # reject
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
