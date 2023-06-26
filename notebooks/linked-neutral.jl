using MultilocusIsland, StatsBase, ColorSchemes, ThreadTools
using Plots, PlotThemes, Distributions, Parameters; theme(:hokusai)
cs = ColorSchemes.viridis

Ls = 1.
L  = 50
s  = Ls/L
h  = 0.5
Ns = 8.
Ne = Ns/s
k  = 5
N  = _Ne2N(Ne, k)
u  = s*0.05
ms = 0.3
A  = Architecture(DipLocus(-s*h, -s), L)
M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), ms*s, ones(L))
P,_= fixedpointit(M, ones(1));
p  = P[end,1,1]
pq = P[end,1,2]

_,Q= simulate(M, 21000, zeros(L), drop=1000, thin=2)
p_ = 1-mean(Q)

# Add a linked neutral locus, and assume it is maintained at p=0.5 on the
# mainland
rs = 5   # strength of recombination relative to selection
A2 = Architecture([DipLocus(0., 0.); A.loci], [rs*s; fill(0.5, L)])
M2 = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A2), ms*s, [0.5; ones(L)])

_, Q2 = simulate(M2, 21000, [rand(); zeros(L)], drop=1000, thin=2)

# Compute expected heterozygosity
P2, _ = MultilocusIsland.fixedpointit_linkage(M2, ones(L+1));
p2    = P2[end,1,1]
pq2   = P2[end,1,2]
EFst  = 1 - P2[end,1,2] / 0.25

p2_ = 1-mean(Q2[:,2:end])
pqs = Q2[:,1] .* (1 .- Q2[:,1])
Fst = 1 .- pqs/0.25
histogram(Fst)

function melinked(locus, p, pq, r, q=1-p)
	@unpack sa, sb, m = locus
	N = (sa*p + sb*pq)
	D = m + r - sa + (2sa + sb)*q - 2sb*pq
	m*(1 + N/D)
end

mes = MultilocusIsland.gff(M, p, pq)
mel = melinked(mes[1], p, pq, rs*s)
EFst2 = 1/(1 + 2Ne*mel)
vline!([EFst], color=:red)
vline!([EFst2], color=:blue)
vline!([mean(Fst)], color=:black)

# Allele frequency distribution
histogram(Q2[:,1], normalize=true)
xs = expectedsfs(M2, p2, pq2)
plot!(xs[1]...)

A3 = Architecture([DipLocus(0., 0.)])
M3 = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A3), ms*s, [0.5])
P3, _ = fixedpointit(M3, ones(1));
p3    = P3[end,1,1]
pq3   = P3[end,1,2]
xs = expectedsfs(M3, p3, pq3)
plot!(xs[1]...)

h=0.5
Xs = tmap(enumerate([0.1, 0.5, 1, 10, 25])) do (i,rs)
    A  = Architecture(DipLocus(-s*h, -s), L)
    A2 = Architecture([DipLocus(0., 0.); A.loci], [rs*s; fill(0.5, L)])
    M2 = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A2), ms*s, [0.5; ones(L)])
    P2, _ = MultilocusIsland.fixedpointit_linkage(M2, ones(L+1));
    p2    = P2[end,:,1]
    pq2   = P2[end,:,2]
    xs2   = MultilocusIsland.expectedsfs_linkage(M2, p2, pq2)
    _, Q2 = simulate(M2, 210000, [rand(); zeros(L)], drop=10000, thin=5)
    rs, xs2, Q2
end

plot()
map(enumerate(Xs)) do (i,(rs,xs2,Q2))
    plot!(xs2[1]..., color=i, label="\$r/s = $rs\$")
    plot!(sfs(Q2[:,1], step=0.05), color=i, marker=true, lw=1, ms=3, alpha=0.5, label="")
end
A3 = Architecture([DipLocus(0., 0.)])
M3 = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A3), ms*s, [0.5])
P3, _ = fixedpointit(M3, ones(1));
p3    = P3[end,1,1]
pq3   = P3[end,1,2]
xs3   = expectedsfs(M3, p3, pq3)
plot!(xs3[1]..., color=:black, lw=2, alpha=0.5, label="single locus")
plot!(xlabel="\$p\$", ylabel="\$\\phi(p)\$")

