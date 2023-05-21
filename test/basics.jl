using MultilocusIsland

L  = 20
Ls = 1.
s  = Ls/L
h  = 0.9
τ  = 0.8
N  = 220
k  = 5
u  = s*0.005
A  = Architecture(HapDipLocus(-s*(1-τ), -s*h*τ, -s*τ), L) 

xs = map(0.05:0.1:0.85) do ms
    M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), ms*s, ones(L)) 
    _, P = simulate(M, 21000, zeros(L), drop=1000, thin=2)
    Q, _ = fixedpointit(M, [1.0])
    mean(P), 1-Q[end,1,1]
end

