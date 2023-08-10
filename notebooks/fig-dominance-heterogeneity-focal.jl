using MultilocusIsland, StatsBase, ColorSchemes, ThreadTools, Printf, QuadGK
using Plots, PlotThemes, Distributions, Parameters; theme(:hokusai)
using LogExpFunctions
cs = ColorSchemes.viridis

L  = 100
Ls = 1.5 
κ  = 1
s̄  = Ls/L
λ  = κ/s̄
Ns = 20
Ne = Ns/s̄
k  = 5
N  = _Ne2N(Ne, k)
u  = s̄ * 0.005
sd = Gamma(κ, 1/λ)
dfe1 = IndependentDFE(sd, Beta(1, 1))

# Contrast weakly and strongly selected recessives and dominants
ss = quantile(sd, [0.05, 0.25, 0.5, 0.75, 0.95])
hs = [0.0, 0.5, 1.0]
mss1 = 0.05:0.05:1.25
mss2 = 0.05:0.01:1.25
nrep = 20
results = map(hs) do h
    map(ss) do s
        @info s, h
        ys1 = map(mss1) do ms
            # random architectures
            map(1:nrep) do _
                A  = Architecture([[randlocus(dfe1) for i=1:L-1] ; DipLocus(-s*h, -s)])
                M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), ms*s̄)
                P,_= fixedpointit(M, ones(L));
                pm = P[end,end,1]
            end
        end
        # average architecture
        ys2 = map(mss2) do ms
            A  = Architecture([[DipLocus(-s̄/2, -s̄) for i=1:L-1] ; DipLocus(-s*h, -s)])
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), ms*s̄)
            P, _ = fixedpointit(M, ones(2))
            pm = P[end,end,1]
        end
        (h, s, ys1, ys2)
    end |> x->vcat(x...)
end



P = plot()
for (i,xss) in enumerate([xs1,xs2,xs3,xs4])
    map(zip(xss,mss)) do (xs,ms)
        scatter!(repeat([ms], length(xs)), xs, alpha=0.3, ms=2, color=i)
    end
    plot!(mss, map(mean, xss), color=i)
end
hline!([0.62], color=:black, ls=:dot)
plot(P)
plot!(mss2, xs5, color=:black)
plot!(mss2, xs6, color=:gray)
