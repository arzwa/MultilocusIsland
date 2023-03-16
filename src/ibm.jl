# Individual-based simulations
Genome{T}     = Vector{T} 
Population{T} = Vector{Genome{T}}

function simulate(model, n; pinit=zeros(nloci(model)), drop=0, thin=1)
    simulate(default_rng(), model, n, pinit, drop, thin)
end

function simulate(rng, model, n, pinit, drop=0, thin=1)
    pop = initialize_ibm(rng, model, pinit)
    ps  = allelefreqs(pop.island)
    Ps  = Matrix{eltype(ps)}(undef, n, nloci(model))
    Ps[1,:] = ps
    for i=2:n
        pop = generation(rng, model, pop)
        Ps[i,:] = allelefreqs(pop.island)
    end
    return pop, Ps[drop+1:thin:end,:]
end

struct MainlandIslandPop{M,I}
    mainland::M
    island::Population{I}
end

function initialize_ibm(rng::AbstractRNG, M::HapDipMainlandIsland, pinit)
    mainland = mainlandpop(M.y)     
    island   = map(_->map(i->rand(rng) < pinit[i] ? 1 : 0, 1:nloci(M)), 1:M.N)
    return MainlandIslandPop(mainland, island)
end

# mainland populations
function mainlandpop(y)
    all(y .== 1) && (return FixedMainland(ones( Int, length(y))))
    all(y .== 0) && (return FixedMainland(zeros(Int, length(y))))
    return PolymorphicMainland(y)
end

struct FixedMainland{T}
    genotype::T
end

struct PolymorphicMainland{T}
    p::T
end

function sample(rng::AbstractRNG, model::FixedMainland, n) 
    [copy(model.genotype) for i=1:n]   # XXX copy!
end

function sample(rng::AbstractRNG, model::PolymorphicMainland, n)
    map(_->[rand(rng) < p ? 0 : 1 for p in model.p], 1:n)
end

function haploidfitness(A::Architecture{T,V}, pop::Population) where {T,V}
    w = Vector{T}(undef, length(pop))
    for i=1:length(pop)
        w[i] = haploidfitness(A, pop[i])
    end
    return lognormalize(w)
end

# Haplodiplontic generation with mainland island migration
# XXX could be faster with some preallocation I guess...
generation(model, pop) = generation(default_rng(), model, pop)
function generation(rng::AbstractRNG, 
        model::HapDipMainlandIsland, 
        pop  ::MainlandIslandPop{M,I}) where {M,I}
    @unpack N, k, m, arch = model
    Nk = N*k
    # migration from mainland
    island = migration(rng, pop, m)
    # calculate haploid fitness
    wg = haploidfitness(arch, island)
    # sample gametes according to haploid fitness
    idxa = sample(rng, 1:N, Weights(wg), 2Nk, replace=true)
    dips = Vector{Tuple{Genome{I},Genome{I}}}(undef, Nk)
    ws = Vector{Float64}(undef, Nk)
    for i=1:Nk
        p1 = island[idxa[i]] 
        p2 = island[idxa[Nk+i]]
        gt = p1 .+ p2  # diploid genotype
        ws[i] = diploidfitness(arch, gt)
        dips[i] = (p1, p2)
    end
    # sample meiospores according to diploid fitness
    ws   = lognormalize(ws)
    idxb = sample(rng, 1:Nk, Weights(ws), N, replace=true)
    island′ = Vector{Genome{I}}(undef, N)
    for (i,j) in enumerate(idxb)
        haplotypes = dips[j]
        island′[i] = meiosis(rng, arch, haplotypes)
    end
    mutation!(rng, island′, model)
    return reconstruct(pop, island=island′)
end

#function generation(rng::AbstractRNG, model::HapMainlandIsland, pop::MainlandIslandPop)
#    @unpack N, m, arch = model
#    # migration from mainland
#    island = migration(rng, pop, m)
#    # calculate haploid fitness
#    wg = lognormalize(map(x->haploidfitness(arch, x), island))
#    # sample next generation according to haploid fitness
#    idxa = sample(rng, 1:N, Weights(wg), 2N, replace=true)
#    off = map(1:N) do i
#        meiosis(rng, island[idxa[i]] .+ island[idxa[N+i]])
#    end
#    mutation!(rng, off, model)
#    return reconstruct(model, island=island)
#end

function migration(rng::AbstractRNG, metapop::MainlandIslandPop, m)
    @unpack mainland, island = metapop
    N = length(island)
    nmigrants = min(N, rand(rng, Poisson(m * N)))  # `min` is neede with strong migration!
    migrants  = sample(rng, mainland, nmigrants)
    residents = copy.(sample(rng, island, N - nmigrants, replace=false)) 
    # XXX copy! replace=false!
    return vcat(migrants, residents)
end

function mutation!(rng::AbstractRNG, pop, model::MainlandIslandModel)
    @unpack N, u = model
    L = length(pop[1])
    n = rand(rng, Poisson(u*L*N))
    ix = sample(rng, 1:N, n, replace=true)
    jx = sample(rng, 1:L, n, replace=true)
    for (i,j) in zip(ix, jx)
        x = copy(pop[i])
        pop[i][j] = (pop[i][j] + 1) % 2
    end
end

function meiosis(rng::AbstractRNG, arch, hs)
    k = rand(rng, 1:2)  # parent to start from
    L = length(hs[1])
    g = similar(hs[1]) 
    for i=1:L
        g[i] = hs[k][i]
        k = rand(rng) < arch.rrate[i] ? 1 + k % 2 : k
    end
    return g
end

allelefreqs(pop, ploidy=1) = reduce(.+, pop) ./ (ploidy*length(pop)) 

