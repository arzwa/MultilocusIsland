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

sasb(l::HapDipLocus) = (sa=l.s1 + l.s01, sb=l.s11 - 2l.s01)

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

# haplodiplontic mean fitness ∑ᵢpᵢeˢⁱ(∑ⱼpⱼeˢⁱʲ), reduces properly to diploid
# and haploid cases when the relevant selection coeffs are 0.
meanfitness(l::HapDipLocus, p, q=1-p) = 
    p*(p + exp(l.s01)*q) + exp(l.s1)*q*(exp(l.s01)*p + exp(l.s11)*q)

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
    R    ::Matrix{V}  # recombination rates between all loci
    function Architecture(
            loci::Vector{HapDipLocus{T}}, 
            rrate::Vector{V}) where {T,V}
        new{T,V}(loci, rrate, rrates(rrate))
    end
end
Architecture(loci) = Architecture(loci, fill(0.5, length(loci)))
Architecture(locus::HapDipLocus, L::Int) = Architecture([locus for i=1:L], fill(0.5, L))
Architecture(locus::HapDipLocus, L::Int, r) = Architecture([locus for i=1:L], fill(r, L))

Base.length(A::Architecture) = length(A.loci)
Base.show(io::IO, A::Architecture) = write(io, "Architecture[L=$(length(A)) loci]")
Base.getindex(A::Architecture, i::UnitRange)  = Architecture(A.loci[i], A.rrate[i]) 
Base.getindex(A::Architecture, i::Int) = A.loci[i]

function Base.vcat(A1::Architecture, A2::Architecture)
    Architecture(vcat(A1.loci, A2.loci), vcat(A1.rrate, A2.rrate))
end

haploidfitness(A::Architecture, g) = haploidfitness(A.loci, g)
diploidfitness(A::Architecture, g) = diploidfitness(A.loci, g)

function summarize(A::Architecture)
    L = length(A)
    w = countmap(A.loci)
    I = indexmap(A.loci)
    l = unique(A.loci)
    γ = [w[k] for k in l] 
    K = length(l)
    (loci=l, w=w, γ=γ, K=K, I=I, L=L) 
end
