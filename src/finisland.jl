# The finite island-model, and individual-based simulation thereof.
# - Migration is assumed to work with contribution and redistribution from a
#   pool, each deme is therefore characterized by a single parameter, i.e. its
#   contribution to the pool. Directional migration according to a graph is
#   rather tricky to implement...
# - Selection is encoded by the genetic architecture in each local deme. 
# - We assume a biallelic loci and symmetric and equal mutation rates.

abstract type Deme end

# local model
@with_kw struct HapDipDeme{U<:Architecture,T} <: Deme
    N ::Int
    k ::Int
    Ne::T = harmonicmean(N, 2N*k)
    u ::T    # mutation rate
    A ::U    # genetic architecture
end

nloci(d::HapDipDeme) = length(d.A)

# XXX we should assert that the architectures are of the same length...
struct FiniteIslandModel{T<:Deme,V}
    D::Vector{T}
    M::Vector{V}
end

Base.getindex(M::FiniteIslandModel, i) = M.D[i]
ndeme(M::FiniteIslandModel) = length(M.D)
nloci(M::FiniteIslandModel) = nloci(M[1])

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

function allelefreqs(pops::Vector{<:Population})
    cat(allelefreqs.(pops)..., dims=2)
end

function initialize_ibm(rng::AbstractRNG, M::FiniteIslandModel, p0)
    map(k->initialize_island(rng, M[k].N, nloci(M[k]), p0[k]), 1:ndeme(M))
end

function initialize_island(rng::AbstractRNG, N, L, p0)
    map(_->map(i->rand(rng) < p0[i] ? 1 : 0, 1:L), 1:N)
end

generation(model, pop) = generation(default_rng(), model, pop)
function generation(
        rng  ::AbstractRNG, 
        model::FiniteIslandModel, 
        pops ::Vector{Population{T}}) where T
    D = ndeme(model)
    # migration
    pops = migration(rng, pops, model.M)
    # within deme updates
    pops_ = similar.(pops)
    for k=1:D
        pops_[k] = generation(rng, model[k], pops[k])
    end
    return pops_
end

# This should be shared with mainland-island model.
function generation(
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

# Every deme contributes a Poisson number of individuals to a common pool,
# which are then redistributed to replenish the local populations.
function migration(rng::AbstractRNG, pops::Vector{T}, M) where T
    D  = length(pops)
    Ns = length.(pops)
    Ms = [min(Ns[i], rand(rng, Poisson(Ns[i]*M[i]))) for i=1:D]
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




