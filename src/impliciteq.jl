# Numerical integration to obtain expected allele frequencies, assuming the
# implicit equation in 𝔼[p].
# We use the integration by parts (IBP) approach suggested by Himani.
# For multiclass barriers, we need to solve a system of nonlinear equations,
# each of which involves a numerical integration.
# We also need to compute 𝔼[pq] (i.e. expected heterozygosity) under Wright's
# distribution.

# The gene flow factor at a single locus of a particular class due to the other
# loci of the same and other classes.
function _gff(xs, w, L, j)
    g = 0.
    for i=1:length(xs)
        g += i == j ? (w[i]-1) * xs[i] : w[i] * xs[i]
    end
    return exp(2g)
end

# Gather all parameters for the single-locus Wright distribution models
# assocaietd with each class of loci. This also involves computing the
# effective migration rates for all classes of loci.
function classparams(M, classes, p)
    @unpack m, u = M
    @unpack m, loci, K, L, γ, y = classes
    xs = map(zip(p, loci, γ, y)) do (pj, lj, wj, yj)
        @unpack s1, s01, s11 = lj
        aj = s1 + s01 + (s11 - 2s01)*(1-pj) 
        (yj-(1-pj))*aj
    end
    map(1:K) do j
        @unpack s1, s01, s11 = loci[j]
        g = _gff(xs, γ, L, j)
        (sa=s1+s01, sb=s11-2s01, N=getNe(M), m=m*g, u=u, pm=1-y[j])
    end
end

# with heterozygosity
function classparams(M, classes, p, pq)
    @unpack m, u = M
    @unpack m, loci, K, L, γ, y = classes
    xs = map(zip(p, pq, loci, γ, y)) do (pj, pqj, lj, wj, yj)
        @unpack s1, s01, s11 = lj
        # exp(2L[sa(qₘ - Eq) + sb(E[pq] - pₘEq)])
        sa = s1 + s01 
        sb = s11 - 2s01
        qj = 1-pj
        sa*(yj - qj) + sb*(pqj - (1-yj)*qj) 
    end
    map(1:K) do j
        @unpack s1, s01, s11 = loci[j]
        g = _gff(xs, γ, L, j)
        (sa=s1+s01, sb=s11-2s01, N=getNe(M), m=m*g, u=u, pm=1-y[j])
    end
end

# sum of square residuals
ssr(x) = sum(x .^ 2)

# Compute the expected SFS for each class of loci conditioning on the gff being
# at its expected value? I think this is what's in Himani's fig 1C is?
function expectedsfs(M, q, pq; step=0.05, f=identity, kwargs...)
    classes = summarize_arch(M)
    θ = classparams(M, classes, q, pq)
    x = (step/2):step:(1-step/2)
    map(1:length(θ)) do j
        Z = Zfun(θ[j]; kwargs...)  # normalizing constant
        y = map(x->Zxfun(θ[j], x, x+step; kwargs...)/(Z*step), 0:step:1-step)
        x, f.(y)
    end
end

expectedsfs(M, Q::Array{T,3}; kwargs...) where T =
    expectedsfs(M, Q[end,:,1], Q[end,:,2]; kwargs...)

# with E[pq]
function fixedpointit(M::MainlandIslandModel, p::AbstractVector; tol=1e-9, kwargs...)
    criterion = false
    classes = summarize_arch(M)
    θ   = classparams(M, classes, p)
    pq  = [Epqfun(θ[j]; kwargs...) for j=1:length(θ)]
    ps  = [p]
    pqs = [pq]
    while !criterion
        θ = classparams(M, classes, ps[end], pqs[end])
        Ep  = similar(p)
        Epq = similar(pq)
        for j=1:length(θ)
            Ep[j]  = Epfun( θ[j]; kwargs...)
            Epq[j] = Epqfun(θ[j]; kwargs...)
        end
        ϵ = Ep .- ps[end]
        push!(ps, Ep)
        push!(pqs, Epq)
        criterion = ssr(ϵ) < tol 
    end
    _ps  = permutedims(hcat(ps...))
    _pqs = permutedims(hcat(pqs...))
    return cat(_ps, _pqs, dims=3), θ
end

# computes the gffs for each unique locus
function gff(M, p, pq)
    classes = summarize_arch(M)
    θ = classparams(M, classes, p, pq)
end

# Fixed point iteration with an arbitrary linkage map
function _gff(A::Architecture, p, pq, y, j)
    @unpack loci, rrate = A
    rs = rrates(rrate, j)
    x = 0.
    for i=1:length(p)
        i == j && continue
        @unpack s1, s01, s11 = loci[i]
        sa = s1 + s01
        sb = s11 - 2s01
        qi = 1 - p[i]
        x += (sa*(y[i] - qi) + sb*(pq[i] - (1-y[i])*qi)) / rs[i] 
    end
    return exp(x)  # XXX factor 2 disappears here with linkage -> in the rᵢ
end

function locusparams(M, p, pq)
    @unpack m, u, y = M
    map(1:length(M.arch)) do i
        @unpack s1, s01, s11 = M.arch[i]
        g = _gff(M.arch, p, pq, y, i)
        (sa=s1+s01, sb=s11-2s01, N=getNe(M), m=m*g, u=u, pm=1-y[i])
    end
end

function fixedpointit_linkage(M::MainlandIslandModel, p::AbstractVector; tol=1e-9, kwargs...)
    criterion = false
    θ   = locusparams(M, p, p .* (1 .- p))
    pq  = [Epqfun(θ[j]; kwargs...) for j=1:length(θ)]
    ps  = [p]
    pqs = [pq]
    while !criterion
        θ   = locusparams(M, ps[end], pqs[end])
        Ep  = similar(p)
        Epq = similar(pq)
        for j=1:length(θ)
            Ep[j]  = Epfun( θ[j]; kwargs...)
            Epq[j] = Epqfun(θ[j]; kwargs...)
        end
        ϵ = Ep .- ps[end]
        push!(ps, Ep)
        push!(pqs, Epq)
        criterion = ssr(ϵ) < tol 
    end
    _ps  = permutedims(hcat(ps...))
    _pqs = permutedims(hcat(pqs...))
    return cat(_ps, _pqs, dims=3)
end



# The exponent of p and (1-p) resp. in Wright's distribution, +1
Afun(N, u, m, pm) = 2N*(u + m*pm)
Bfun(N, u, m, pm) = 2N*(u + m*(1-pm))

# The selection part of Wright's distribution
Cfun(N, p, sa, sb) = exp(-N*p*(2sa + sb*(2-p)))

# The single-locus (log)probability density function
ϕfun(p, N, u, m, pm, sa, sb) = 
    p^(Afun(N,u,m,pm)-1) * (1-p)^(Bfun(N,u,m,pm)-1) * Cfun(N,p,sa,sb)
lϕfun(p, N, u, m, pm, sa, sb) = 
    (Afun(N,u,m,pm)-1)*log(p) + (Bfun(N,u,m,pm)-1)*log(1-p) - N*p*(2sa+sb*(2-p))

function singlelocus_sfs(N, u, m, pm, sa, sb; Δ=0.02, kwargs...)
    Z = log(Zfun((N=N, u=u, m=m, pm=pm, sa=sa, sb=sb); kwargs...))
    ps = (Δ/2):Δ:(1-Δ/2)
    xs = map(p->exp(lϕfun(p, N, u, m, pm, sa, sb) - Z), ps)
    ps, reverse(xs ./ sum(xs))
end

function fpfun(p, θ)
	@unpack sa, sb, N, m, u, pm = θ  # these are the parameters we need
    B = Bfun(N, u, m, pm)
    C = Cfun(N, p, sa, sb)
    -(1-p)^(B-1)*(N * (-p*sb + 2sa + sb*(2-p)) + (B-1)/(1-p)) * C 
end

function gpfun(p, θ)
	@unpack sa, sb, N, m, u, pm = θ
    A = Afun(N, u, m, pm)
    C = Cfun(N, p, sa, sb)
    p^(A-1)*(N*(p*sb - 2sa + sb*(p-2)) + (A-1)/p) * C 
end

function gspfun(p, θ)
	@unpack sa, sb, N, m, u, pm = θ
    A = Afun(N, u, m, pm)
    C = Cfun(N, p, sa, sb)
    p^A*(N*(p*sb - 2sa + sb*(p-2)) + A/p) * C
end

# Normalizing constant of Wright's distribution, using IBP
# could just use Z0fun/Z1fun/Zxfun below (but implemented this one first...)
function Zfun(θ; kwargs...)
	@unpack sa, sb, N, m, u, pm = θ
    A = Afun(N, u, m, pm)
    B = Bfun(N, u, m, pm)
	c = 0.5^(A + B - 1) * exp(-N*(sa + 3sb/4)) * (1/A + 1/B)
	d, _ = quadgk(p->((    p^A)/A)*fpfun(p, θ), 0, 0.5; kwargs...)
	e, _ = quadgk(p->(((1-p)^B)/B)*gpfun(p, θ), 0.5, 1; kwargs...)
	c - d + e
end

# 0 to x integral
function Z0fun(θ, x; kwargs...)
	@unpack sa, sb, N, m, u, pm = θ
    A = Afun(N, u, m, pm)
    B = Bfun(N, u, m, pm)
    a = (1/A)*(x)^A*(1-x)^(B-1)*Cfun(N, x, sa, sb)
	b, _ = quadgk(p->((p^A)/A)*fpfun(p, θ), 0, x; kwargs...)
    return a - b
end

# x to 1 integral
function Z1fun(θ, x; kwargs...)
	@unpack sa, sb, N, m, u, pm = θ
    A = Afun(N, u, m, pm)
    B = Bfun(N, u, m, pm)
    a = (1/B)*(1-x)^B*x^(A-1)*Cfun(N, x, sa, sb)
    b, _ = quadgk(p->(((1-p)^B)/B)*gpfun(p, θ), x, 1; kwargs...)
    return a - b
end

# General x to y integral
function Zxfun(θ, x, y; kwargs...)
	@unpack sa, sb, N, m, u, pm = θ
    x ≈ 0 && y ≈ 1 && return Zfun(θ; kwargs...)
    x ≈ 0 && return Z0fun(θ, y; kwargs...)
    y ≈ 1 && return Z1fun(θ, x; kwargs...)
    I, _ = quadgk(p->ϕfun(p, N, u, m, pm, sa, sb), x, y; kwargs...)
    return I
end

function Yfun(θ; kwargs...)
	@unpack sa, sb, N, m, u, pm = θ
    A = Afun(N, u, m, pm)
    B = Bfun(N, u, m, pm)
	C(p) = Cfun(N, p, sa, sb)
	L, _ = quadgk(p->p^A * (1-p)^(B-1) * C(p), 0., 0.5; kwargs...)
	R1 = (1/B)*(1/2)^(B+A)*C(0.5)
	R2, _ = quadgk(p->(((1-p)^B)/B)*gspfun(p, θ), 0.5, 1; kwargs...)
	return L + R1 + R2
end

function Epqfun(θ; kwargs...)
	@unpack sa, sb, N, m, u, pm = θ
    A = Afun(N,u,m,pm)
    B = Bfun(N,u,m,pm)
    N, _ = quadgk(p -> p^A * (1-p)^B * Cfun(N,p,sa,sb), 0., 1.)
    Z = Zfun(θ; kwargs...)
    return N/Z
end

# Expected value
Epfun(θ; kwargs...) = Yfun(θ; kwargs...) / Zfun(θ; kwargs...)

