# Individual-based simulations
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
    island::I
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

# Haplodiplontic generation with mainland island migration
generation(model, pop) = generation(default_rng(), model, pop)
function generation(rng::AbstractRNG, model::HapDipMainlandIsland, pop::MainlandIslandPop)
    @unpack N, k, m, arch = model
    Nk = N*k
    # migration from mainland
    island = migration(rng, pop, m)
    # calculate haploid fitness
    wg = map(x->haploidfitness(arch, x), island) |> lognormalize
    # sample gametes according to haploid fitness
    idxa = sample(rng, 1:N, Weights(wg), 2Nk, replace=true)
    dips = map(1:Nk) do i
        off = island[idxa[i]] .+ island[idxa[Nk+i]] 
        off, diploidfitness(arch, off)
    end
    # sample meiospores according to diploid fitness
    ws   = last.(dips) |> lognormalize
    idxb = sample(rng, 1:Nk, Weights(ws), N, replace=true)
    smpl = first.(dips)[idxb]
    island′ = map(x->meiosis(rng, x), smpl)
    mutation!(rng, island′, model)
    return reconstruct(pop, island=island′)
end

function generation(rng::AbstractRNG, model::HapMainlandIsland, pop::MainlandIslandPop)
    @unpack N, m, arch = model
    # migration from mainland
    island = migration(rng, pop, m)
    # calculate haploid fitness
    wg = lognormalize(map(x->haploidfitness(arch, x), island))
    # sample next generation according to haploid fitness
    idxa = sample(rng, 1:N, Weights(wg), 2N, replace=true)
    off = map(1:N) do i
        meiosis(rng, island[idxa[i]] .+ island[idxa[N+i]])
    end
    mutation!(rng, off, model)
    return reconstruct(model, island=island)
end

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

meiosis(rng::AbstractRNG, diploid) = map(l->segregate(rng, l), diploid)

function segregate(rng::AbstractRNG, l::Int)
    l == 0 ? 0 : (l == 2 ? 1 : (rand(rng) < 0.5 ? 0 : 1))
end

allelefreqs(pop, ploidy=1) = reduce(.+, pop) ./ (ploidy*length(pop)) 

