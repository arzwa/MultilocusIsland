using MultilocusIsland, StatsBase, ColorSchemes
using Plots, PlotThemes, Distributions, Parameters; theme(:hokusai)
cs = ColorSchemes.viridis

Ls = 1.3
L  = 100
s  = Ls/L
h  = 0.9
Ns = 8.
Ne = Ns/s
k  = 5
N  = _Ne2N(Ne, k)
u  = s*0.005
ms = 0.3
A  = Architecture(DipLocus(-s*h, -s), L)
M  = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A, y=ones(L))
P,_= fixedpointit(M, ones(1));
p  = P[end,1,1]
pq = P[end,1,2]

AA = deepcopy(A)
push!(AA, DipLocus(0., 0.), NaN)
AA.rrate[end-1] = 0.02s
M2  = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=AA, y=[ones(L); 0.5])
_, Q = simulate(M2, 11000)

pqs = Q[:,end] .* (1 .- Q[:,end])
Fst = 1 .- pqs/0.25
histogram!(Fst)

function me(M, p, pq)
	classes = summarize_arch(M)
	θ = MultilocusIsland.classparams(M, classes, p, pq)
end

function melinked(locus, p, pq, r, q=1-p)
	@unpack sa, sb, m = locus
	N = (sa*p + sb*pq)
	D = m + r - sa + (2sa + sb)*q - 2sb*pq
	m*(1 + N/D)
end

mes = me(M, p, pq)

mel = melinked(mes[1], p, pq, 5s)

1/(1 + 2Ne*mel)
