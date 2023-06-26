# The finite island-model, and individual-based simulation thereof.
# - Migration is assumed to work with contribution and redistribution from a
#   pool, each deme is therefore characterized by a single parameter, i.e. its
#   contribution to the pool. Directional migration according to a graph is
#   rather tricky to implement...
# - Selection is encoded by the genetic architecture in each local deme. 
# - We assume biallelic loci and symmetric and equal mutation rates.

# Data structures for IBM
Genome{T}     = Vector{T} 
Population{T} = Vector{Genome{T}}

# Local deme model
# ================
abstract type Deme end

"""
    HapDipDeme

Model for a haplodiplontic deme, with a discrete and synchronized alternation
of generations between N haploids and Nk diploids.
"""
@with_kw struct HapDipDeme{U<:Architecture,T} <: Deme
    N ::Int
    k ::Int
    Ne::T = harmonicmean(N, 2N*k)
    u ::T  # mutation rate
    A ::U  # genetic architecture
end

nloci(d::HapDipDeme) = length(d.A)


# Finite islands model
# ====================
"""
    FiniteIslandModel
"""
struct FiniteIslandModel{T<:Deme,V}
    D::Vector{T}
    m::Vector{V}  # migration rates
    function FiniteIslandModel(d::Vector{T}, m::Vector{V}) where {T,V}
        @assert length(unique([length(di.A) for di in d])) == 1
        @assert length(d) == length(m)
        return new{T,V}(d, m)
    end
end

Base.getindex(M::FiniteIslandModel, i) = M.D[i]
ndeme(M::FiniteIslandModel) = length(M.D)
nloci(M::FiniteIslandModel) = nloci(M[1])


# Mainland-island type model
# ==========================
# Note, we use a dedicated type for the mainland population, this enables more
# efficient simulation (no random variables have to be generated in the
# important special cases of a mainland fixed for either allele!)
abstract type Mainland end

struct FixedMainland{T,V} <: Mainland
    genotype::T
    p::Vector{V}
end

struct PolymorphicMainland{T} <: Mainland
    p::Vector{T}
end

sample(rng::AbstractRNG, M::FixedMainland, n) = [copy(M.genotype) for i=1:n]
sample(rng::AbstractRNG, M::PolymorphicMainland, n) = 
    map(_->[rand(rng) < p ? 1 : 0 for p in M.p], 1:n)

"""
    MainlandIslandModel 
"""
struct MainlandIslandModel{DD<:Deme,MM<:Mainland}
    D::DD
    m::Float64     # migration rate
    mainland::MM
end

nloci(M::MainlandIslandModel) = nloci(M.D)

function MainlandIslandModel(D, m, y::Vector{<:Real})
    if all(y .== 1) 
        M = FixedMainland(y, ones(Int, length(y)))
    elseif all(y .== 0) 
        M = FixedMainland(y, zeros(Int, length(y)))
    else
        M = PolymorphicMainland(y)
    end
    return MainlandIslandModel(D, m, M)
end


# Simulation for the finite islands model
function simulate(model::FiniteIslandModel, n, p0; drop=0, thin=1)
    simulate(default_rng(), model, n, p0, drop, thin)
end

function simulate(rng, model::FiniteIslandModel, n, p0, drop=0, thin=1)
    pop = initialize_ibm(rng, model, p0)
    Ps   = Array{Float64,3}(undef, n, nloci(model), ndeme(model)) 
    Ps[1,:,:] = allelefreqs(pop)
    @showprogress 1 "[Running simulation]" for i=2:n
        pop = generation(rng, model, pop)
        Ps[i,:,:] = allelefreqs(pop)
    end
    return pop, Ps[drop+1:thin:end,:,:]
end

# simulation for the mainland island model
function simulate(model::MainlandIslandModel, n, p0; drop=0, thin=1)
    simulate(default_rng(), model, n, p0, drop, thin)
end

function simulate(rng, model::MainlandIslandModel, n, p0, drop=0, thin=1)
    pop = initialize_island(rng, model.D.N, length(model.D.A), p0)
    Ps   = Array{Float64,2}(undef, n, nloci(model)) 
    Ps[1,:] = allelefreqs(pop)
    @showprogress 1 "[Running simulation]" for i=2:n
        pop = generation(rng, model, pop)
        Ps[i,:] = allelefreqs(pop)
    end
    return pop, Ps[drop+1:thin:end,:]
end

allelefreqs(pops::Vector{<:Population}) = cat(allelefreqs.(pops)..., dims=2)

initialize_ibm(rng::AbstractRNG, M, p0) = p0

function initialize_ibm(rng::AbstractRNG, M::FiniteIslandModel, p0::Matrix{Float64})
    map(k->initialize_island(rng, M[k].N, nloci(M[k]), p0[k,:]), 1:ndeme(M))
end

function initialize_island(rng::AbstractRNG, N, L, p0)
    map(_->map(i->rand(rng) < p0[i] ? 1 : 0, 1:L), 1:N)
end

generation(model, pop) = generation(default_rng(), model, pop)

# Finite islands model
function generation(
        rng  ::AbstractRNG, 
        model::FiniteIslandModel, 
        pops ::Vector{Population{T}}) where T
    D = ndeme(model)
    # migration
    pops = migration(rng, pops, model)
    # within deme updates
    pops_ = similar.(pops)
    for k=1:D
        pops_[k] = reproduction(rng, model[k], pops[k])
    end
    return pops_
end

# Mainland-island model
function generation(
        rng  ::AbstractRNG, 
        model::MainlandIslandModel, 
        pop  ::Population{T}) where T
    # migration
    pop = migration(rng, pop, model)
    return reproduction(rng, model.D, pop)
end

# This should be shared with mainland-island model.
function reproduction(
        rng  ::AbstractRNG, 
        model::HapDipDeme, 
        pop  ::Population{T}) where T
    @unpack N, k, A = model
    Nk = N*k
    # calculate haploid fitness
    wg = haploidfitness(A, pop)
    # sample gametes according to haploid fitness
    idxa = sample(rng, 1:N, Weights(wg), 2Nk, replace=true)
    dips = Vector{Tuple{Genome{T},Genome{T}}}(undef, Nk)
    ws   = Vector{Float64}(undef, Nk)
    for i=1:Nk
        p1 = pop[idxa[i]] 
        p2 = pop[idxa[Nk+i]]
        gt = p1 .+ p2  # diploid genotype
        ws[i] = diploidfitness(A, gt)
        dips[i] = (p1, p2)
    end
    # sample meiospores according to diploid fitness
    ws   = lognormalize(ws)
    idxb = sample(rng, 1:Nk, Weights(ws), N, replace=true)
    pop_ = Vector{Genome{T}}(undef, N)
    for (i,j) in enumerate(idxb)
        haplotypes = dips[j]
        pop_[i] = meiosis(rng, A, haplotypes)
    end
    mutation!(rng, pop_, model)
    return pop_
end

function haploidfitness(A::Architecture{T,V}, pop::Population) where {T,V}
    w = Vector{T}(undef, length(pop))
    for i=1:length(pop)
        w[i] = haploidfitness(A, pop[i])
    end
    return lognormalize(w)
end

# Every deme contributes a Poisson number of individuals to a common pool,
# which are then redistributed to replenish the local populations.
function migration(rng::AbstractRNG, pops::Vector{T}, M::FiniteIslandModel) where T
    D  = length(pops)
    Ns = length.(pops)
    Ms = [min(Ns[i], rand(rng, Poisson(Ns[i]*M.m[i]))) for i=1:D]
    newpop = similar.(pops)
    shuffle!.(rng, pops)
    pool = vcat([copy(pops[i][1:Ms[i]]) for i=1:D]...)::T
    shuffle!(rng, pool)
    k = 1
    for i=1:D
        newpop[i][1:Ms[i]] .= pool[k:k+Ms[i]-1]   # random from migrant pool
        newpop[i][Ms[i]+1:end] .= copy(pops[i][Ms[i]+1:end])  # stayers
        k += Ms[i]
    end
    return newpop
end

function migration(rng::AbstractRNG, pop, M::MainlandIslandModel)
    N = length(pop)
    nmigrants = min(N, rand(rng, Poisson(M.m * N)))  # `min` is needed with strong migration!
    migrants  = sample(rng, M.mainland, nmigrants)
    residents = copy.(sample(rng, pop, N - nmigrants, replace=false)) 
    # XXX copy! replace=false!
    return vcat(migrants, residents)
end

function mutation!(rng::AbstractRNG, pop, model)
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

allelefreqs(pop::Population{T}, ploidy=1) where T = 
    reduce(.+, pop) ./ (ploidy*length(pop)) 

