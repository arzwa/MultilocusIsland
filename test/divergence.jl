using MultilocusIsland, Random, Plots, PlotThemes, StatsBase, ThreadTools; theme(:hokusai)

# divergent selection, with dominance 
L  = 25
Ls = 1.2
s  = Ls/L
N  = 50
k  = 5
h1 = 0.50
h2 = 0.50
D1 = 1
D2 = 1
A1 = Architecture(DipLocus(-s*h1, -s), L)
A2 = Architecture(DipLocus( s*h2,  s), L)
As = [[A1 for i=1:D1] ; [A2 for i=1:D2]]
d  = D1 + D2
m  = fill(0.2, d) .* s
D  = [MultilocusIsland.HapDipDeme(N=N, k=k, A=A, u=0.001*s) for A in As]
MM = MultilocusIsland.FiniteIslandModel(D, m)
p0 = permutedims(hcat(zeros(L, D1) .+ 1e-3, ones(L, D2) .- 1e-3))

X, P1 = simulate(MM, 110000, p0, drop=10000, thin=10)
div = sum(abs.(P1[:,:,1] .- P1[:,:,2]), dims=2)

# expand the architecture

A1b = vcat(A1, Architecture(DipLocus(   0., 0.), 3L))
A2b = vcat(A2, Architecture(DipLocus( s*h2,  s), 3L))

init = deepcopy(X)
for pop in init
    for i=1:length(pop)
        pop[i] = [pop[i] ; zeros(Int, 3L)]
    end
end
        
Asb = [[A1b for i=1:D1] ; [A2b for i=1:D2]]
Db  = [MultilocusIsland.HapDipDeme(N=N, k=k, A=A, u=0.001*s) for A in Asb]
MMb = MultilocusIsland.FiniteIslandModel(Db, m)
_, P2 = simulate(MMb, 110000, init, drop=0, thin=1)

divb = sum(abs.(P2[:,:,1] .- P2[:,:,2]), dims=2)

plot([div; divb])


