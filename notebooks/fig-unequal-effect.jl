
using ThreadTools, Printf, Plots, PlotThemes, Serialization, MCMCChains
using Base.Iterators, Parameters, Random, StatsBase, MultilocusIsland
using ColorSchemes, Distributions, DataFrames, CSV
theme(:hokusai)

# We want to sample appropriate dominance coefficients for a given τ. There are
# obviously multiple ways to do this.

# This function assumes a single shared selection coefficient, and a branching
# process approximation for the fixation probability. Under these assumptions
# the distribution of the fixed `h` does not depend on `s` nor `Ne`. It does
# give a pronounced Haldane's sieve effect (much more so than when we use a
# diffusion-based fixation probability for not very large populations).
# Note that the distribution has a closed form, but it is more straightforward
# to sample from a sufficiently fine discretization.
function sample_h1(hdist, τ, n)
    ps = map(h->(1-τ*(1-h))*pdf(hdist, h), 0:0.01:1)
    ps ./= sum(ps)
    1 .- rand(DiscreteNonParametric(0:0.01:1, ps), n)
end

# **Note** that the dominance coefficients `h` for beneficials fixed on the
# island in allopatry correspond to dominance coefficients `1-h` of invading
# (ancestral) alleles from the mainland in the mainland-island model.

# Let us start by contrasting a haploid and diploid life cycle (`τ=0` vs.
# 'τ=1`). We do the computation with variable `h` and a shared `h` set to the
# expected value.

Ne2N(Ne, k) = ceil(Int, Ne/(2k) + Ne)
k  = 5
t  = 0.5
Ns = 4.  
s  = 0.02
N  = Ne2N(Ns/s, k) 
Ls = 1.5
L  = ceil(Int, Ls/s)
u  = 0.005*s
a, b = 1, 1
hd = Beta(a,b)  # distribution of `h` for new mutations
sd = Exponential(s)  # distribution of `h` for new mutations
hs = sample_h1(hd, t, L)
ss = rand(sd, L)
arch = Architecture([HapDipLocus(-s*(1-t), -s*h*t, -s*t) for (s,h) in zip(ss,hs)])
GS(L) = GibbsSampler([UnitIntervalProposal() for i=1:L])

Xs = map(0:0.05:1) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
    Q = fixedpointit(M, zeros(L))[end,:,1]
end

plot(1 .- permutedims(hcat(Xs...)))

M = HapDipMainlandIsland(N=N, k=k, m=0.5*s, u=u, arch=arch)
_, Q1 = simulate(M, 210000, drop=10000, thin=100)
Q2 = gibbs(M, GS(L), rand(L), 55000, drop=5000)
Q3 = fixedpointit(M, zeros(L))[end,:,1]
qq = expectedsfs(M, Q3, step=0.05)

PP = plot(title="\$L=$L, L\\bar{s}=$Ls, N_e\\bar{s}=$Ns, \\tau=$t, \\alpha=$a, \\beta=$b\$")
map(enumerate(1:15:L)) do (i,j)
    plot!(sfs(Q2[:,j], step=0.05, f=log10), color=i, label=string(arch[j]))
    x, y = qq[j]
    plot!(x, log10.(y), color=i, lw=6, alpha=0.2, label="")
    scatter!(sfs(Q1[:,j], step=0.05, f=log10), color=i, markerstrokecolor=i, marker=true, label="")
end
plot(PP, xlabel="\$q\$", ylabel="\$\\log_{10}\\phi\$", size=(500,250))

plot(sort(Q3))
plot!(sort(vec(mean(Q2, dims=1))))
plot!(sort(vec(mean(Q1, dims=1))))
