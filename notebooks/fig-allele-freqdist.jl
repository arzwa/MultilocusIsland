using ThreadTools, Printf, Plots, PlotThemes, Serialization, MCMCChains
using Base.Iterators, Parameters, Random, StatsBase, MultilocusIsland
using ColorSchemes, Distributions, DataFrames, CSV
theme(:hokusai)

Ne2N(Ne, k) = ceil(Int, Ne/(2k) + Ne)

s  = 0.02
ms = [0.162, 0.486, 0.742, 0.813]
Ns = [2, 4, 8, 16]
k  = 5
N  = Ne2N.(Ns ./ s, Ref(k)) 
t  = 0.0
L  = 40
u  = 0.005*s
h  = 0.5

P1 = plot()
Xs = map(zip(1:4, ms, N)) do (i, m, N)
    A = Architecture(HapDipLocus(-s*(1-t), -s*h*t, -s*t), L)
    G = GibbsSampler(L)
    M = HapDipMainlandIsland(N=N, k=k, m=m*s, u=u, arch=A)
    q = fixedpointit(M, [0.])[end,1,1]
    Q = gibbs(M, GibbsSampler(L), rand(L), 5500, drop=500)
    x, y = expectedsfs(M, [q], step=0.05)[1]
    plot!(x, log10.(y), color=i, alpha=0.2, lw=6)
    x, y = sfs(vec(Q), step=0.05, f=log10)
    plot!(x, y, color=i, ls=:dot)
    _, Q = simulate(M, 110000, drop=10000, thin=10)
    x, y = sfs(vec(Q), step=0.05, f=log10)
    plot!(x, y, color=i, markerstrokecolor=i, marker=true, ylim=(-6,0))
end
plot(P1)

savefig("/home/arthur_z/vimwiki/build/img/2023-03-30/himani1ci-b.svg")
    

P2 = plot()
Ys = map(zip(1:4, ms, N)) do (i, m, N)
    A = Architecture(HapDipLocus(-s*(1-t), -s*h*t, -s*t), L)
    G = GibbsSampler(L)
    M = HapDipMainlandIsland(N=N, k=k, m=m*s, u=u, arch=A)
    q = fixedpointit(M, [0.])[end,1,1]
    Q = gibbs(M, GibbsSampler(L), rand(L), 5500, drop=500)
    x, y = expectedsfs(M, [q], step=0.05)[1]
    plot!(x, log10.(y), color=i, alpha=0.2, lw=6)
    x, y = sfs(vec(Q), step=0.05, f=log10)
    plot!(x, y, color=i, ls=:dot)
    Ne = Ns[i] ./ s
    M = HapMainlandIsland(N=Ne, m=m*s, u=u, arch=A)
    _, Q = simulate(M, 110000, drop=10000, thin=10)
    x, y = sfs(vec(Q), step=0.05, f=log10)
    plot!(x, y, color=i, markerstrokecolor=i, marker=true, ylim=(-6,0))
end
plot(P2)


A = Architecture(HapDipLocus(-s*(1-t), -s*h*t, -s*t), L)
G = GibbsSampler(L)
M = HapDipMainlandIsland(N=N[1], k=k, m=ms[1]*s, u=u, arch=A)
q = fixedpointit(M, [0.])[end,1,1]

Q2 = gibbs(M, GibbsSampler(L), rand(L), 5500, drop=500)

_, Q = simulate(M, 110000, drop=10000, thin=10) 

x, y = expectedsfs(M, [q], step=0.05, f=log10)[1]
plot(x, y)
x2, y2 = sfs(vec(Q2), step=0.05, f=log10)
plot!(x2, y2)
x, y = sfs(vec(Q), step=0.05, f=log10)
plot!(x, y)
