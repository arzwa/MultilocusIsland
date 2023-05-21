using MultilocusIsland, Random, Plots, PlotThemes, StatsBase, ThreadTools; theme(:hokusai)

s = 0.02
L = 40
d = 3
A = [Architecture([HapLocus(-s) for i=1:L], fill(0.5, L)) for i=1:d]
D = [MultilocusIsland.HapDipDeme(N=200, k=10, A=Architecture(HapLocus(-s), L), u=0.005*s) for i=1:d]
M = [s*0.1, s*0.1, s*0.1]
model = MultilocusIsland.FiniteIslandModel(D, M)
p0 = [rand(L) ./ 4 for i=1:L]
_, P = simulate(model, 10000, p0)


# divergent selection, with dominance 
L  = 25
Ls = 1
s  = Ls/L
N  = 50
k  = 5
h1 = 0.10
h2 = 0.10
D1 = 1
D2 = 1
A1 = Architecture(DipLocus(-s*h1, -s), L)
A2 = Architecture(DipLocus( s*h2,  s), L)
As = [[A1 for i=1:D1] ; [A2 for i=1:D2]]
d  = D1 + D2
m  = fill(0.2, d) .* s
D  = [MultilocusIsland.HapDipDeme(N=N, k=k, A=A, u=0.005*s) for A in As]
MM = MultilocusIsland.FiniteIslandModel(D, m)
p0 = permutedims(hcat(zeros(L, D1) .+ 1e-3, ones(L, D2) .- 1e-3))

_, P = simulate(MM, 110000, p0, drop=10000, thin=10)
QQ, HH = fixedpointit(MM, p0)

mP1 = mean(P[:,:,1], dims=2)
mP2 = mean(P[:,:,2], dims=2)
plot(mP1)
plot!(mP2)
hline!(1 .- QQ[:,1,end])

xx = expectedsfs(MM, QQ[:,:,end], HH[:,:,end], f=log10)

plot(xx[1][1]); plot!(xx[2][1])
scatter!(sfs(vec(P[:,:,1]), f=log10))
scatter!(sfs(vec(P[:,:,2]), f=log10))

MultilocusIsland.ϕ(MM, rand(d,L), d, 1)
G = MultilocusIsland.GibbsSampler(MM)
PP = MultilocusIsland.gibbs(MM, G, p0, 10000)

P1 = plot(sfs(vec(PP[:,1,:]), f=log10), 
          title="\$Ls=$Ls, L=$L, D_1 = $D1, D_2=$D2, h=$h1, N=$N, m/s=$(m[1]/s)\$")
plot!(sfs(vec(PP[:,d,:]), f=log10))
scatter!(sfs(vec(P[:,:,1]), f=log10))
scatter!(sfs(vec(P[:,:,d]), f=log10))


# as the number of demes grows, it should start working?
L  = 25
Ls = 1
s  = Ls/L
N  = 50
k  = 5
h1 = 0.10
h2 = 0.10
D1 = 25
D2 = 25
A1 = Architecture(DipLocus(-s*h1, -s), L)
A2 = Architecture(DipLocus( s*h2,  s), L)
As = [[A1 for i=1:D1] ; [A2 for i=1:D2]]
d  = D1 + D2
m  = fill(0.2, d) .* s
D  = [MultilocusIsland.HapDipDeme(N=N, k=k, A=A, u=0.005*s) for A in As]
MM = MultilocusIsland.FiniteIslandModel(D, m)

# Gibbs? 
MultilocusIsland.ϕ(MM, rand(d,L), d, 1)
G = MultilocusIsland.GibbsSampler(MM)
p0 = permutedims(hcat(zeros(L, D1) .+ 1e-3, ones(L, D2) .- 1e-3))
PP = MultilocusIsland.gibbs(MM, G, p0, 10000)
_, P = simulate(MM, 21000, p0, drop=1000)

P2 = plot(sfs(vec(PP[:,1,:]), f=log10), 
          title="\$Ls=$Ls, L=$L, D_1 = $D1, D_2=$D2, h=$h1, N=$N, m/s=$(m[1]/s)\$")
plot!(sfs(vec(PP[:,d,:]), f=log10))
scatter!(sfs(vec(P[:,:,1]), f=log10))
scatter!(sfs(vec(P[:,:,d]), f=log10))


plot(P1, P2, legend=false, size=(800,250), ms=2, titlefont=8)


ps = mean(P, dims=2)

plot(ps[:,:,1], ylim=(0,1))
plot!(ps[:,:,2])

D = abs.(P[:,:,1] .- P[:,:,2])

plot()
for i=1:10:L
    plot!(D[1:1:end,i])
end
plot!()


# Swamping thresholds in the symmetric model
L  = 30
N  = 220
k  = 5
s  = 0.05
hs = [0.5, 0.75, 1.0, 1.5]
mss= 0.05:0.2:1.95
Ys = map(hs) do h
    @info h
    ys = map(mss) do ms
        m  = fill(ms, 2) .* s
        A1 = Architecture(DipLocus(-s*h, -s), L)
        A2 = Architecture(DipLocus( s*h,  s), L)
        push!(A1, DipLocus(0., 0.), 0.5)
        push!(A2, DipLocus(0., 0.), 0.5)
        p0 = [zeros(L+1), ones(L+1)]
        D  = [MultilocusIsland.HapDipDeme(N=N, k=k, A=A, u=0.005*s) for A in [A1,A2]]
        MM = MultilocusIsland.FiniteIslandModel(D, m)
        _, P = simulate(MM, 11000, p0, drop=1000)
        (ms, P)
    end
    (h, ys)
end

P1 = plot(xlabel="\$m/s\$", ylabel="\$\\tilde{p}\$", legend=:bottomleft)
P2 = plot(xlabel="\$m/s\$", ylabel="\$\\Delta\$", legend=:topright)
P3 = plot(xlabel="\$m/s\$", ylabel="\$F_{ST}\$", legend=:bottomright)
for (i, (h, X)) in enumerate(Ys)
    xs = first.(X)
    ys = map(P->mean(P[:,1:L,2]), last.(X))
    plot!(P1, xs, ys, label="\$h = $h\$", color=i, lw=2, marker=true, markerstrokecolor=i)
    ys = map(P->mean(P[:,1:L,1]), last.(X))
    plot!(P1, xs, 1 .- ys, label="", lw=2, alpha=0.3, color=i, marker=true, markerstrokecolor=i)
    ys = map(last.(X)) do P
        D = P[:,1:L,2] .- P[:,1:L,1]
        mean(D)
    end
    plot!(P2, xs, ys, label="\$h=$h\$", color=i, lw=2, marker=true, markerstrokecolor=i)
    ys = map(last.(X)) do P
        p1 = P[:,end,1]
        p2 = P[:,end,2]
        p̄1 = mean(p1)
        p̄2 = mean(p2)
        p̄ = mean([p1 ; p2])
        Fst = 1 - 0.5*(p̄1*(1-p̄1) + p̄2*(1-p̄2))/p̄*(1-p̄)
    end
    plot!(P3, xs, ys, label="\$h=$h\$", color=i, lw=2, marker=true, markerstrokecolor=i)
end
plot(P1, P2, P3, ylim=(0,1), layout=(1,3), size=(800,200))

## Infinitish island model
L  = 30
D1 = 45
D2 = 5
N  = 220
k  = 5
s  = 0.05
h  = 1.0
mss= 0.05:0.2:1.25
X  = map(mss) do ms
    m  = fill(ms, D1+D2) .* s
    A1 = Architecture(DipLocus(-s*h, -s), L)
    A2 = Architecture(DipLocus( s*h,  s), L)
    push!(A1, DipLocus(0., 0.), 0.5)
    push!(A2, DipLocus(0., 0.), 0.5)
    D  = [MultilocusIsland.HapDipDeme(N=N, k=k, A=A1, u=0.005*s) for i=1:D1]
    D  = [D ; [MultilocusIsland.HapDipDeme(N=N, k=k, A=A2, u=0.005*s) for i=1:D2]]
    p0 = [[zeros(L+1) for i=1:D1] ; [ones(L+1) for i=1:D2]]
    MM = MultilocusIsland.FiniteIslandModel(D, m)
    _, P = simulate(MM, 11000, p0, drop=1000)
    (ms, P)
end

P1 = plot(xlabel="\$m/s\$", ylabel="\$\\tilde{p}\$", legend=:bottomleft)
P2 = plot(xlabel="\$m/s\$", ylabel="\$\\Delta\$", legend=:topright)
P3 = plot(xlabel="\$m/s\$", ylabel="\$F_{ST}\$", legend=:bottomright)
xs = first.(X)
ys = map(P->mean(P[:,1:L,1:D1]), last.(X))
plot!(P1, xs, 1 .- ys)
ys = map(P->mean(P[:,1:L,D1+1:end]), last.(X))
plot!(P1, xs, ys, label="")
ys = map(last.(X)) do P
    D = abs(mean(P[:,1:L,1:D1]) - mean(P[:,1:L,D1+1:end]))
end
plot!(P2, xs, ys)
ys = map(last.(X)) do P
    p1 = P[:,end,1]
    p2 = P[:,end,2]
    p̄1 = mean(p1)
    p̄2 = mean(p2)
    p̄ = mean([p1 ; p2])
    Fst = 1 - 0.5*(p̄1*(1-p̄1) + p̄2*(1-p̄2))/p̄*(1-p̄)
end
plot!(P3, xs, ys)
plot(P1, P2, P3, ylim=(0,1), layout=(1,3), size=(800,200))




