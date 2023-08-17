# # MultilocusIsland.jl

# This julia package implements:
#
# 1. Individual-based simulations of a multilocus, biallelic, mainland-island
#    population genetics model for a haplodiplontic/haploid/diploid life cycle.
# 2. Numerical tools to calculate effective migration rates and approximate
#    allele frequency distributions in the above model using diffusion theory.
#
# To install the package:
#
# 1. Install `julia` from [https://julialang.org/](https://julialang.org/).
# 2. Open up a terminal and start a REPL by entering `julia` in the prompt.
# 3. Type `]` to enter the julia package manager.
# 4. Type `add https://github.com/arzwa/MultilocusIsland` to install the package.
#
# After executing these steps, one should be able to load the package using
# `using MultilocusIsland` in the julia REPL.

# ## Example

using MultilocusIsland, Distributions, Random; Random.seed!(12)

# This defines a genetic architecture with `L` biallelic loci with
# mean selection coefficient `s̄`, dominance coefficients distributed
# according to a Beta(2,2) distribution, and relative strength of haploid vs.
# diploid selection drawn from a uniform distribution:
function random_locus(s̄)
    s = rand(Exponential(s̄))
    h = rand(Beta(2,2))
    t = rand()
    return HapDipLocus(-s*(1-t), -s*h*t, -s*t)
end

L = 50
s̄ = 0.02
A = Architecture([random_locus(s̄) for i=1:L])

# This defines a mainland-island model with `N` haploid individuals per
# generation, `N*k` diploid individuals per generation, mutation rate
# `u`, genetic architecture `A`, haploid migration rate `m1`,
# diploid migration rate `m2` and mainland allele frequencies `y` 
N  = 500
k  = 5
u  = s̄*0.002
m1 = s̄*0.4
m2 = 0.0
y  = ones(L)
M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m1, m2, y)

# Obtain the equilibrium allele frequencies and heterozygosities
P, _ = fixedpointit(M, ones(L))
P
# The second argument to `fixedpointit` is the initial value (locally
# beneficial allele frequencies on the island) for the fixed-point iteration
# algorithm, we use `[1,1,...,1]`, assuming this corresponds to the approach to
# equilibrium starting from secondary contact.
# **Note:** when dealing with homogeneous architectures (i.e. some proportion
# of loci in `A` have identical effects), then the length of the initial vector
# should be the number of *unique* loci in the barrier.

# `P` is a 3D array, `P[:,:,1]` gives the result of the fixed point iteration
# for the allele frequencies, `P[:,:,2]` for the heterozygosities. The columns
# correspond to the different loci, the rows correspond to different iterates
# of the algorithm. The equilibrium allele frequencies hence are:
p = P[end,:,1]

# and the heterozygosities
pq = P[end,:,2]

# We can conduct an individual-based simulation for the same model
genomes, Q = simulate(M, 10000, zeros(L));

# Note that this returns the frequency of the locally deleterious alleles
# (`q`), not locally beneficial as in the above. We can compare the simulations
# against the numerics:
p_ibm = 1 .- vec(mean(Q[5000:5:end,:], dims=1))
collect(zip(p_ibm, p))

# Of course, one should do more/longer simulations.
