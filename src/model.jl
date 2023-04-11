# Population and locus model
# Locus definition
"""
    HapDipLocus(s1, s01, s11)

Haploid-diploid selected locus. Relative fitness of `0` allele is assumed to be
1 in both haploids and diploid homozygotes.
"""
struct HapDipLocus{T}
    s1  :: T  # haploid derived allele selection coefficient
    s01 :: T  # diploid heterozygous effect
    s11 :: T  # diploid homozygous effect
end

Base.show(io::IO, l::HapDipLocus) = write(io, "HapDipLocus$(string(l))")
Base.string(l::HapDipLocus) = @sprintf "(%.3f, %.3f, %.3f)" l.s1 l.s01 l.s11

# We are dealing with fitness on the log-scale. Note that we assume fitness is
# multiplicative across loci.
haploidfitness(l::HapDipLocus, g) = g == 0 ? 0. : l.s1
haploidfitness(l::AbstractVector, g::AbstractVector) = 
    mapreduce(x->haploidfitness(x...), +, zip(l, g))

diploidfitness(l::HapDipLocus, g) = g == 0 ? 0. : g == 1 ? l.s01 : l.s11
diploidfitness(l::AbstractVector, g::AbstractVector) = 
    mapreduce(x->diploidfitness(x...), +, zip(l, g))

# aliases
HapLocus(s) = HapDipLocus(s, 0., 0.)
DipLocus(s01, s11) = HapDipLocus(0., s01, s11)

# Genetic architecture
"""
    Architecture

Genetic architecture, consists of a bunch of loci with certain effects, and a
linkage map (recombination rates). We assume `rrate[i]` gives the **probability
of a crossover** between locus `i` and locus `i+1`. 
"""
struct Architecture{T,V}
    loci ::Vector{HapDipLocus{T}}
    rrate::Vector{V}  # recombination rates between neighboring loci
end
Architecture(loci) = Architecture(loci, fill(0.5, length(loci)))
Architecture(locus::HapDipLocus, L::Int) = Architecture([locus for i=1:L], fill(0.5, L))

Base.length(A::Architecture) = length(A.loci)
Base.show(io::IO, A::Architecture) = write(io, "Architecture[L=$(length(A)) loci]")
Base.getindex(A::Architecture, i)  = A.loci[i] 

function Base.vcat(A1::Architecture, A2::Architecture)
    Architecture(vcat(A1.loci, A2.loci), vcat(A1.rrate, A2.rrate))
end

function Base.push!(A::Architecture, l, r)
    push!(A.loci, l)
    push!(A.rrate,r)
end

haploidfitness(A::Architecture, g) = haploidfitness(A.loci, g)
diploidfitness(A::Architecture, g) = diploidfitness(A.loci, g)


# Mainland - Island type models
abstract type MainlandIslandModel end

Base.show(io::IO, m::T) where T<:MainlandIslandModel = 
    write(io, "$(T.name.name)$((m.N, m.m, m.μ, length(m.island[1])))")  

harmonicmean(x, y) = 1/(1/x + 1/y)

"""
    HapDipMainlandIsland{A}

Haplodiplontic mainland-island model. Note that migration is assumed to occur
only in the haploid phase (think of spore dispersal).
"""
@with_kw struct HapDipMainlandIsland{A} <: MainlandIslandModel
    N::Int
    k::Int
    m::Float64
    u::Float64
    arch::A  # genetic architecture
    y::Vector{Float64} = ones(length(arch)) # mainland allele frequencies
    Ne::Float64        = harmonicmean(N, 2N*k)
end
# note `y` are the frequencies of the derived alleles, i.e. what we would call
# q_mainland (frequency of the selected alleles on the mainland). 

getNe(m::HapDipMainlandIsland) = m.Ne

# Not retested
# Haploid only -- no `k`
"""
    HapMainlandIsland{A}

Haplontic mainland-island model.
"""
@with_kw struct HapMainlandIsland{A} <: MainlandIslandModel
    N::Int
    m::Float64
    u::Float64
    arch::A      # genetic architecture
    y::Vector{Float64} = ones(length(arch))
end

getNe(m::HapMainlandIsland) = m.N

nloci(M::MainlandIslandModel) = length(M.arch)

function sfs(ps; step=0.02, f=identity)
    hs = fit(Histogram, ps, 0:step:1+step) 
    hs = StatsBase.normalize(hs, mode=:probability)
    es = collect(hs.edges)[1][2:end-1] .- step/2
    ws = hs.weights[1:end-1]    # make sure we add the fixed states to the last bin 
    ws[end] += hs.weights[end] 
    es, f.(ws)
end

# For these things and the expected values, it would be nicer of the
# architecture was maintained in a more informative type, which we specialize
# for the different applications (mainly IBM/Gibbs vs. numerics)
function summarize_arch(M::MainlandIslandModel)
    # get the different classes of loci
    L = length(M.arch)
    w = countmap(M.arch.loci)
    I = indexmap(M.arch.loci)
    l = unique(M.arch.loci)
    γ = [w[k] for k in l] 
    y = [M.y[I[k]] for k in l] 
    K = length(l)
    return (m=M.m, loci=l, L=L, γ=γ, K=K, y=y)
end
