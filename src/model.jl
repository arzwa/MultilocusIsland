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

Base.show(io::IO, l::HapDipLocus) = 
    write(io, "Locus(s₁=$(l.s1), s₀₁=$(l.s01), s₁₁=$(l.s11))")

# We are dealing with fitness on the log-scale. Note that we assume fitness is
# multiplicative across loci.
haploidfitness(l::HapDipLocus, g) = g == 0 ? 0. : l.s1
haploidfitness(l::AbstractVector, g::AbstractVector) = 
    mapreduce(x->haploidfitness(x...), +, zip(l, g))

diploidfitness(l::HapDipLocus, g) = g == 0 ? 0. : g == 1 ? l.s01 : l.s11
diploidfitness(l::AbstractVector, g::AbstractVector) = 
    mapreduce(x->diploidfitness(x...), +, zip(l, g))

# Mainland - Island type models
abstract type MainlandIslandModel end

Base.show(io::IO, m::T) where T<:MainlandIslandModel = 
    write(io, "$(T.name.name)$((m.N, m.m, m.μ, length(m.island[1])))")  

harmonicmean(x, y) = 1/(1/x + 1/y)

"""
    HapDipMainlandIsland{A}

Haplodiplontic mainland-island model.
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

# Not retested
# Haploid only -- no `k`
"""
    HapMainlandIsland{A}

Haplontic mainland-island model.
"""
@with_kw struct HapMainlandIsland{A,I,M} <: MainlandIslandModel
    N::Int
    m::Float64
    u::Float64
    arch::A      # genetic architecture
    y::Vector{Float64} = ones(length(arch))
end

nloci(M::MainlandIslandModel) = length(M.arch)

function sfs(ps; step=0.02, f=identity)
    hs = fit(Histogram, ps, 0:step:1) 
    hs = StatsBase.normalize(hs, mode=:probability)
    es = collect(hs.edges)[1][2:end] .- step/2
    ws = f.(hs.weights)
    es, ws
end

# For these things and the expected values, it would be nicer of the
# architecture was maintained in a more informative type, which we specialize
# for the different applications (mainly IBM/Gibbs vs. numerics)
function summarize_arch(M::MainlandIslandModel)
    # get the different classes of loci
    L = length(M.arch)
    w = countmap(M.arch)
    I = indexmap(M.arch)
    l = unique(M.arch)
    γ = [w[k] for k in l] 
    y = [M.y[I[k]] for k in l] 
    K = length(l)
    return (m=M.m, loci=l, L=L, γ=γ, K=K, y=y)
end
