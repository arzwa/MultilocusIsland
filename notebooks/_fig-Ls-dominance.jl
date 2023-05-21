# Compare approximations to individual-based sims
using MultilocusIsland, StatsBase, ColorSchemes, ThreadTools
using Plots, PlotThemes, Distributions, Parameters; theme(:hokusai)
cs = ColorSchemes.viridis

# 1. Ls appreciable, decreasing L, increasing s (within reasonable bounds) for 
# various levels of dominance. This is the main thing...

# Say Ls = 1, Nes = 5
Ls  = 1.
Nes = 5.
k   = 5
s   = [0.2, 0.1, 0.05, 0.01]

Ne  = Nes ./ s
NN  = Ne2N.(Ne, k)
LL  = ceil.(Int, Ls ./ s)
hs  = [0., 0.5, 1.]

xs  = vec(collect(Iterators.product(zip(NN, LL), hs)))
ms1 = 0:0.01:0.6
ms2 = 0.05:0.1:0.55

X1 = tmap(xs) do ((N, L), h)
    @info (N,L,h)
    s = Ls/L
    A = Architecture(DipLocus(-s*h, -s), L)
    y1 = map(ms2) do ms
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=s*0.005, arch=A)
        _, P = simulate(M, 110000, drop=10000, thin=10) 
        ms, mean(P)
    end
    N, L, h, y1
end

X2 = map(xs) do ((N, L), h)
    @info (N,L,h)
    s = Ls/L
    A = Architecture(DipLocus(-s*h, -s), L)
    y2 = map(ms1) do ms
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=s*0.005, arch=A)
        P,_= fixedpointit(M, [1.])
        pm = P[end,1,1]
        ms, pm
    end
    N, L, h, y2
end

map(1:length(X1)) do i
    (N,L,h,y1) = X1[i]
    xs, ys = first.(y1), 1 .- last.(y1)
    scatter(xs, ys, title="\$N=$N, L=$L, h=$h\$", legend=false, markerstrokecolor=1,
            ms=4, color=:black, ylim=(0,1), xlim=(0,0.61))
    (N,L,h,y2) = X2[i]
    xs, ys = first.(y2), last.(y2)
    plot!(xs, ys, color=:gray, lw=3, alpha=0.5)
end |> x->plot(x..., size=(750,420), titlefont=8, xlabel="\$m/s\$", ylabel="\$\\tilde{p}\$")

savefig("notebooks/img/Ls-dominance.svg")


# show the bistability 
Ls  = 1.
Nes = 10.
k   = 5
L   = 50
s   = Ls/L
u   = s*0.01
h   = 1.
Ne  = Nes/ s
N   = _Ne2N(Ne, k)
A = Architecture(DipLocus(-s*h, -s), L)

ms1 = 0:0.01:0.6
y1 = map(ms1) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A)
    P1,_ = fixedpointit(M, [0.])
    P2,_ = fixedpointit(M, [1.])
    ms, P1[end,1,1], P2[end,1,1]
end

ms3 = 0.1:0.05:0.5
y3 = tmap(ms3) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A)
    n = 0.18 ≤ ms ≤ 0.42 ? 50000 : 20000
    _,P1 = simulate(M, n, pinit=zeros(L), drop=n÷2, thin=2)
    _,P2 = simulate(M, n, pinit=ones(L), drop=n÷2, thin=5)
    @info ms
    ms, mean(P1), mean(P2)
end

ms2 = 0.2:0.1:0.4
y2 = map(ms2) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A)
    PP = plot()
    ps = map(0:0.1:1) do p0
        P1,_ = fixedpointit(M, [p0])
        P1[:,1,1]
    end
    ll = maximum(length.(ps))
    P = plot(framestyle=:default)
    map(ps) do p1
        plot!([p1 ; [p1[end] for i=1:(ll-length(p1))]],
             color=:black, alpha=0.6, legend=false)
    end
    _, xmx = xlims(P)
    _, ymx = ylims(P)
    annotate!(xmx, ymx*1.1, text("\$m/s=$ms\$", 7, :right))
    P
end
xlabel!(y2[end], "iteration")
P2 = plot(y2..., layout=(length(y2),1), size=(160,700), ylabel="\$p_0\$", ytickfont=5, xtickfont=5)
P1 = plot(first.(y1), getindex.(y1,3), label="\$p_0 = 1\$", lw=2, xlabel="\$m/s\$", ylabel="\$\\mathbb{E}[p]\$", 
         title="\$h=1, L=$L, Ls = $Ls, N_es = $Nes\$")
plot!(first.(y1), getindex.(y1,2), label="\$p_0 = 0\$", legend=:topright, lw=2)
scatter!(first.(y3), 1 .- getindex.(y3,2), color=1, markerstrokecolor=1, label="")
scatter!(first.(y3) .+ 0.005, 1 .- getindex.(y3,3), color=2, markerstrokecolor=2, label="")
plot(P1, P2, size=(550,300), layout=grid(1,2,widths=[0.8,0.2]))


