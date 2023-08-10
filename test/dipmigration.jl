using MultilocusIsland, Random, ThreadTools, StatsBase, Plots, PlotThemes; theme(:hokusai)

L  = 40
Ls = 0.8
s  = Ls/L
h  = 0.5
τ  = 0.5
k  = 2
Ns = 8.
Ne = Ns/s
N  = _Ne2N(Ne, k)
u  = s*0.005
A  = Architecture(HapDipLocus(-s*(1-τ), -s*h*τ, -s*τ), L) 

m = 0.35*s
M1    = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m/2, m/2, ones(L)) 
Q1, _ = fixedpointit(M1, [1.0])
_, P1 = simulate(M1, 21000, zeros(L), drop=1000, thin=2)
M2    = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m, 0., ones(L)) 
Q2, _ = fixedpointit(M2, [1.0])
_, P2 = simulate(M2, 21000, zeros(L), drop=1000, thin=2)
M3    = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), 0., m, ones(L)) 
Q3, _ = fixedpointit(M3, [1.0])
X, P3 = simulate(M3, 21000, zeros(L), drop=1000, thin=2)

mean(P1), 1-Q1[end,1,1]
mean(P2), 1-Q2[end,1,1]
mean(P3), 1-Q3[end,1,1]

rng = Random.Xoshiro(123)
pat, mat = MultilocusIsland.haploidphase(rng, M3.D, X)
XX = MultilocusIsland.migration(rng, pat, mat, M3)


ms1 = [0.1, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7]
ms2 = 0.01:0.01:0.8
as  = [0., 0.5, 1]
ys = map(as) do a
    @info a
    ys1 = tmap(ms1) do ms
        M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), a*ms*s, (1-a)*ms*s, ones(L)) 
        _, P = simulate(M, 25000, zeros(L), drop=5000, thin=2)
        1 - mean(P)
    end
    ys2 = map(ms2) do ms
        M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), a*ms*s, (1-a)*ms*s, ones(L)) 
        Q, _ = fixedpointit(M, [1.0])
        Q[end,1,1]
    end
    ys1, ys2
end

plot(title="\$Ls = $Ls, N_es = $Ns, ks = $s, N=$N, k=$k, \\tau=$τ, h=$h\$", titlefont=7)
map(enumerate(zip(as, ys))) do (i,(a,y))
    plot!(ms2, y[2], color=i, label="\$m_1 = $a m, m_2 = $(1-a) m\$")
    scatter!(ms1, y[1], color=i, label="")
end
plot!(legend=:bottomleft, xlabel="\$m/s\$", ylabel="\$\\mathbb{E}[p]\$", size=(350,250))

M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), 0., 0.55*s, ones(L)) 
_, P = simulate(M, 55000, zeros(L), drop=5000, thin=2)

scatter!([0.55], [1 - mean(P[1:1000,:])])


