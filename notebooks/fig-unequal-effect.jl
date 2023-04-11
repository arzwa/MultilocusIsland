
using ThreadTools, Printf, Plots, PlotThemes, Serialization, MCMCChains
using Base.Iterators, Parameters, Random, StatsBase, MultilocusIsland
using ColorSchemes, Distributions, DataFrames, CSV
theme(:hokusai)
cs = ColorSchemes.viridis

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
    ps = map(h->(1-τ*(1-h))*pdf(hdist, h), 0:0.0001:1)
    ps ./= sum(ps)
    1 .- rand(DiscreteNonParametric(0:0.0001:1, ps), n)
end

# **Note** that the dominance coefficients `h` for beneficials fixed on the
# island in allopatry correspond to dominance coefficients `1-h` of invading
# (ancestral) alleles from the mainland in the mainland-island model.
# These will generally be skewed towards 0, i.e. local adaptation being
# dominant (invaders recessive), i.e. the Haldane's sieve effect.

# Let us start by contrasting a haploid and diploid life cycle (`τ=0` vs.
# 'τ=1`). We do the computation with variable `h` and a shared `h` set to the
# expected value.

Ne2N(Ne, k) = ceil(Int, Ne/(2k) + Ne)
k  = 5
t  = 1.0
Ns = 8.  
s  = 0.01
N  = Ne2N(Ns/s, k) 
Ls = 1.0
L  = ceil(Int, Ls/s)
u  = 0.005*s
a, b = 1, 1
hd = Beta(a,b)  # distribution of `h` for new mutations
#sd = Exponential(s)  # distribution of `h` for new mutations
sd = Dirac(s)  # distribution of `h` for new mutations
hs = sample_h1(hd, t, L)
ss = rand(sd, L)
A1 = Architecture([HapDipLocus(-s*(1-t), -s*h*t, -s*t) for (s,h) in zip(ss,hs)])
A2 = Architecture(HapDipLocus(-s*(1-t), -s*mean(hs)*t, -s*t), L)
A3 = Architecture(HapDipLocus(-s*(1-t), -s*t, -s*t), L)
A4 = Architecture(HapDipLocus(-s*(1-t), 0., -s*t), L)

mss = 0:0.02:1
Xs = map(mss) do ms
    M1 = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A1)
    M2 = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A2)
    M3 = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A3)
    M4 = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A4)
    Q1 = fixedpointit(M1, ones(length(unique(hs))))[end,:,1]
    Q2 = fixedpointit(M2, [1.])[end,1,1]
    Q3 = fixedpointit(M3, [1.])[end,1,1]
    Q4 = fixedpointit(M4, [1.])[end,1,1]
    Q1, Q2, Q3, Q4
end

P1 = plot(mss, permutedims(hcat(first.(Xs)...)),
          color=reshape(get.(Ref(cs), hs), 1, length(hs)),
          alpha=0.2, label="", ylabel="\$\\tilde{p}\$", xlabel="\$m/s\$")
plot!(P1, repeat([1.03], length(0:0.01:1)), 0:0.01:1, color=get.(Ref(cs), 0:0.01:1), lw=8, label="")
plot!(P1, mss, getindex.(Xs, 2), color=:white, lw=2, label="")
plot!(P1, mss, getindex.(Xs, 3), color=:black, ls=:dot, label="\$h=1\$")
plot!(P1, mss, getindex.(Xs, 4), color=:black, ls=:dash, label="\$h=0\$", legend=:topright)
P2 = stephist(hs, bins=0:0.1:1, color=:black, legend=false, ylabel="frequency",
              xlabel="\$h\$", normalize=:probability)
plot(P1, P2, size=(500,200))


## Variable s as well
Ne2N(Ne, k) = ceil(Int, Ne/(2k) + Ne)
k  = 5
t  = 1.0
Ns = 8.  
s  = 0.01
N  = Ne2N(Ns/s, k) 
Ls = 1.0
L  = ceil(Int, Ls/s)
u  = 0.005*s
a, b = 1, 1
hd = Beta(a,b)  # distribution of `h` for new mutations
sd = Exponential(s)  # distribution of `h` for new mutations
hs = sample_h1(hd, t, L)
ss = rand(sd, L)
A1 = Architecture([HapDipLocus(-s*(1-t), -s*h*t, -s*t) for (s,h) in zip(ss,hs)])
A2 = Architecture(HapDipLocus(-s*(1-t), mean(-ss .* hs)*t, -s*t), L)
A3 = Architecture(HapDipLocus(-s*(1-t), -s*t, -s*t), L)
A4 = Architecture(HapDipLocus(-s*(1-t), 0., -s*t), L)

sa = mean(-(1-t)*ss .+ (-t * ss .* hs))
sb = mean(-t * ss .- (-2t .* ss .* hs))
S  = -(sb + 2sa) 
H  = 0.5*(1 - sb/S)
A5 = Architecture(HapDipLocus(0., -S*H, -S), L)

mss = 0:0.02:1
Xs = map(mss) do ms
    M1 = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A1)
    M2 = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A2)
    M3 = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A3)
    M4 = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A4)
    M5 = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A5)
    Q1 = fixedpointit(M1, ones(length(unique(hs))))[end,:,1]
    Q2 = fixedpointit(M2, [1.])[end,1,1]
    Q3 = fixedpointit(M3, [1.])[end,1,1]
    Q4 = fixedpointit(M4, [1.])[end,1,1]
    Q5 = fixedpointit(M5, [1.])[end,1,1]
    Q1, Q2, Q3, Q4, Q5
end

Y = permutedims(hcat(first.(Xs)...))
ym = mean(Y, dims=2)
P1 = plot(mss, Y,
          color=reshape(get.(Ref(cs), hs), 1, length(hs)),
          alpha=0.2, label="", ylabel="\$\\tilde{p}\$", xlabel="\$m/s\$")
plot!(P1, mss, ym, color=:red, lw=2, label="")
plot!(P1, repeat([1.03], length(0:0.01:1)), 0:0.01:1, color=get.(Ref(cs), 0:0.01:1), lw=8, label="")
plot!(P1, mss, getindex.(Xs, 2), color=:black, lw=2, label="")
plot!(P1, mss, getindex.(Xs, 3), color=:black, ls=:dot, label="\$h=1\$")
plot!(P1, mss, getindex.(Xs, 4), color=:black, ls=:dash, label="\$h=0\$", legend=:topright)
plot!(P1, mss, getindex.(Xs, 5), color=:blue, label="", legend=:topright)
P2 = stephist(hs, bins=0:0.1:1, color=:black, legend=false, ylabel="frequency",
              xlabel="\$h\$", normalize=:probability)
plot(P1, P2, size=(500,200))

