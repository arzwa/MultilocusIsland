using Random, MultilocusIsland, MCMCDiagnosticTools, StatsBase, ThreadTools
using Plots, PlotThemes; theme(:hokusai)

# sample from a beta distribution using adaptive metropolis
function mhsample(a, b, n)
    d = Beta(a,b)
    prop = AdaptiveProposal(trans = Bijectors.Logit(0,1), stop=10000)
    xs = Vector{Float64}(undef, n)
    xs[1] = rand(d)
    for i=2:n
        _x, q = prop(xs[i-1])
        a = logpdf(d, _x) - logpdf(d, xs[i-1]) + q
        if log(rand()) < a 
            xs[i] = _x 
            accept!(prop)
        else 
            xs[i] = xs[i-1]
        end
    end
    return xs, prop
end

L = 20
s = 0.05
h = 1.5
τ = 1.
N = 200
k = 5
u = s*0.005
mss = 0:0.1:1
mss2 = 0.0:0.02:1.
A = Architecture([HapDipLocus(-s*(1-τ), -s*h*τ, -s*τ) for i=1:L], fill(0.5, L)) 

xs1 = tmap(mss) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, arch=A, u=u)
    G = GibbsSampler([UnitIntervalProposal() for i=1:L])
    _, P = simulate(M, 11000, drop=100, thin=10)
    @info ms
    mean(P)
end

xs2 = tmap(mss) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, arch=A, u=u)
    G = GibbsSampler([UnitIntervalProposal() for i=1:L])
    P, _ = gibbs(M, G, rand(L), 11000, drop=1000)
    @info ms
    mean(P), G
end;

plot(mss, 1 .- first.(xs2))
scatter!(mss, 1 .- xs1)

M = HapDipMainlandIsland(N=N, k=k, m=0.7*s, arch=A, u=u)
G = GibbsSampler([UnitIntervalProposal() for i=1:L]);
ls = map(1:10) do _
    P, l = gibbs(M, G, rand(L), 11000, drop=1000)
    l
end
plot(ls)

PP = ess(reshape(P, (10000,1,20)))
