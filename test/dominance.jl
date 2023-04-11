using MultilocusIsland
using Plots, PlotThemes, Printf, ColorSchemes; theme(:hokusai)
cs = ColorSchemes.viridis

L = 21
s = 0.04
h = 0.0:(1/(L-1)):1
τ = 1.
N = 220
k = 5
u = s*0.005
m = 0.01

# Haplodiplontic model
A = Architecture([HapDipLocus(-s*(1-τ), -s*h[i]*τ, -s*τ) for i=1:L], fill(0.5, L)) 
    
M = HapDipMainlandIsland(N=N, k=k, m=0.01, arch=A, u=u)

ys = map(0:0.01:1) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, arch=A, u=u)
    Eq = fixedpointit(M, ones(L))
    Eq[end,:]
end

plot(0:0.01:1, permutedims(hcat(ys...)), 
     label=reshape(map(x->"\$h=$(@sprintf "%.2f" x)\$", h), 1, length(h)),
     color=reshape(reverse(get.(Ref(cs), h)), 1, length(h)))

_, P = simulate(M, 1000)
G = GibbsSampler([UnitIntervalProposal() for i=1:L])
Q, l = gibbs(M, G, rand(L), 10000)
