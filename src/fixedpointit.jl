# Numerical integration to obtain expected allele frequencies, assuming the
# implicit equation in 𝔼[p].
# We use the integration by parts (IBP) approach suggested by Himani.
# For multiclass barriers, we need to solve a system of nonlinear equations,
# each of which involves a numerical integration.
# We also need to compute 𝔼[pq] (i.e. expected heterozygosity) under Wright's
# distribution.
function getclasses(M::MainlandIslandModel)
    @unpack D, m, mainland = M
    # get the different classes of loci
    θ = summarize(D.A)
    y = [mainland.p[θ.I[k]] for k in θ.loci] 
    return (; θ..., m=M.m, y=y, u=D.u, Ne=D.Ne)
end

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
function classparams(classes, p, pq)
    @unpack m, u, Ne, loci, K, L, γ, y = classes
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
        (sa=s1+s01, sb=s11-2s01, N=Ne, m=m*g, u=u, pm=1-y[j])
    end
end

# sum of square residuals
ssr(x) = sum(x .^ 2)

# Compute the expected SFS for each class of loci conditioning on the gff being
# at its expected value? I think this is what's in Himani's fig 1C is?
function expectedsfs(M, p, pq; step=0.05, f=identity, kwargs...)
    classes = getclasses(M)
    θ = classparams(classes, p, pq)
    x = (step/2):step:(1-step/2)
    map(1:length(θ)) do j
        Z = Zfun(θ[j]; kwargs...)  # normalizing constant
        y = map(x->Zxfun(θ[j], x, x+step; kwargs...)/(Z*step), 0:step:1-step)
        x, f.(y)
    end
end

expectedsfs(M, Q::Array{T,3}; kwargs...) where T =
    expectedsfs(M, Q[end,:,1], Q[end,:,2]; kwargs...)

function fixedpointit(
        M::MainlandIslandModel, 
        p::AbstractVector; 
        tol=1e-9, kwargs...)
    criterion = false
    classes = getclasses(M)
    θ   = classparams(classes, p, p .* (1 .- p))
    pq  = [Epqfun(θ[j]; kwargs...) for j=1:length(θ)]
    ps  = [p]
    pqs = [pq]
    while !criterion
        θ = classparams(classes, ps[end], pqs[end])
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
    classes = getclasses(M)
    θ = classparams(classes, p, pq)
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

function locusparams(M::MainlandIslandModel, p, pq)
    @unpack u, A = M.D
    y = M.mainland.p
    map(1:length(A)) do i
        @unpack s1, s01, s11 = A[i]
        g = _gff(A, p, pq, y, i)
        (sa=s1+s01, sb=s11-2s01, N=M.D.Ne, m=M.m*g, u=u, pm=1-y[i])
    end
end

function fixedpointit_linkage(
        M::MainlandIslandModel, 
        p::AbstractVector; 
        tol=1e-9, kwargs...)
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
    return cat(_ps, _pqs, dims=3), θ
end


# A wild attempt for the finite-island model
function migrantpool(m, p)
    (K,L) = size(p)
    y = zeros(K, L)
    for k=1:K, j=1:K
        y[k,:] .+= (p[k,:] .* m[k])
    end
    return y ./ sum(m)
end

function locusparams(M::FiniteIslandModel, p, pq)
    @unpack m, D = M
    Y = migrantpool(m, p)
    map(1:ndeme(M)) do k
        @unpack A, Ne, u = D[k]
        map(1:length(A)) do i
            @unpack s1, s01, s11 = A[i]
            g = _gff(A, p[k,:], pq[k,:], Y[k,:], i)
            (sa=s1+s01, sb=s11-2s01, N=Ne, m=m[k]*g, u=u, pm=1-Y[k,i])
        end
    end
end

function fixedpointit(
        M::FiniteIslandModel, 
        p::AbstractMatrix; 
        α=0.1, tol=1e-9, maxit=100, kwargs...)
    criterion = false
    K,L = size(p)
    θ   = locusparams(M, p, p .* (1 .- p))
    pq  = [Epqfun(θ[k][j]; kwargs...) for k=1:K, j=1:L]
    ps  = [p]
    pqs = [pq]
    while !criterion
        θ   = locusparams(M, ps[end], pqs[end])
        Ep  = similar(p)
        Epq = similar(pq)
        for k=1:K, j=1:L
            Ep[k,j]  =  ps[end][k,j] * (1-α) + Epfun( θ[k][j]; kwargs...) * α
            Epq[k,j] = pqs[end][k,j] * (1-α) + Epqfun(θ[k][j]; kwargs...) * α
        end
        ϵ = Ep .- ps[end]
        push!(ps, Ep)
        push!(pqs, Epq)
        criterion = ssr(ϵ) < tol || length(ps) > maxit
    end
    cat(ps..., dims=3), cat(pqs..., dims=3)
end

function expectedsfs(
        M::FiniteIslandModel, 
        p::AbstractMatrix, 
        pq::AbstractMatrix; 
        step=0.05, f=identity, 
        kwargs...)
    K,L = size(p)
    θ = locusparams(M, p, pq)
    x = (step/2):step:(1-step/2)
    map(1:K) do k
        map(1:L) do j
            Z = Zfun(θ[k][j]; kwargs...)  # normalizing constant
            y = map(x->Zxfun(θ[k][j], x, x+step; kwargs...)/(Z*step), 0:step:1-step)
            x, f.(y)
        end
    end
end