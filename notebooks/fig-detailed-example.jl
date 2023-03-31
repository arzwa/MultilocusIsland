using ThreadTools, Printf, Plots, PlotThemes, Serialization, MCMCChains
using Parameters, Random, StatsBase, MultilocusIsland
theme(:hokusai)

# It takes quite a while for the IBMs to reach equilibrium it seems. The MRF
# seems to do quite a good job... Should make this more precise...
# There is a strong correlation between different loci.
Ne2N(Ne, k) = ceil(Int, Ne/(2k) + Ne)
k  = 5
h  = 1. 
t  = 1.
Ns = 4.  
s  = 0.02
N  = Ne2N(Ns/s, k) 
Ls = 1.5
L  = ceil(Int, Ls/s)
u  = 0.005*s
ngen = 305000
drop = 5000
thin = 50
arch = Architecture([HapDipLocus(-s*(1-t), -s*h*t, -s*t) for i=1:L])
GS(L) = GibbsSampler([UnitIntervalProposal() for i=1:L])
mmax = 0.4

X1 = tmap(0:0.05:mmax) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
    _, Q1  = simulate(M, ngen, drop=drop, thin=thin) 
    q1 = mean(Q1)
    ms, q1, Q1
end

X2 = tmap(0:0.01:mmax) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
    Q2 = gibbs(M, GS(L), rand(length(arch)), 10500, drop=500)
    q2 = mean(Q2)
    ss = vec(ess(Chains(Q2))[:,2]) 
    ms, q2, Q2, ss
end

X3 = tmap(0:0.001:mmax) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
    q3 = fixedpointit(M, [1.0])[end,:,:]
    q4 = fixedpointit(M, [0.0])[end,:,:]
    ms, q3[1], q4[1], q3[2], q4[2]
end

#serialize("data/h0.5-example.jls", (X1, X2, X3))
serialize("data/h1-example.jls", (X1, X2, X3))

X1, X2, X3 = deserialize("data/h1-example.jls")

P1 = plot(title="\$L=$L, Ls=$Ls, N_es=$Ns, h=$h\$")
plot!(first.(X3), 1 .- getindex.(X3,2), label="E+", alpha=0.5, lw=3)
plot!(first.(X3), 1 .- getindex.(X3,3), label="E-", alpha=0.5, lw=3,
     marker=false,  xlabel="\$m/s\$", ylabel="\$p\$")
plot!(first.(X2), 1 .- getindex.(X2,2), label="MRF", color=:black, alpha=0.4, marker=true, ms=1)
scatter!(first.(X1), markerstrokecolor=:black,
             1 .- map(X->mean(X[1500:end,:]), getindex.(X1,3)), color=:black,
             label="IBM", legend=:bottomleft, ms=3)

    
M = HapDipMainlandIsland(N=N, k=k, m=0.25*s, u=u, arch=arch)
_, Q1  = simulate(M, 110000, pinit=ones(L), drop=10000, thin=10) 
_, Q2  = simulate(M, 110000, pinit=rand(L), drop=10000, thin=10) 
_, Q3  = simulate(M, 110000, pinit=zeros(L), drop=10000, thin=10) 

mean(1 .- Q1)
mean(1 .- Q2)
mean(1 .- Q3)

QQ1 = 1 .- Q1 
P2 = plot(QQ1[:,1], xlabel="generation / 10", ylabel="\$p\$",
          title="\$m/s=$(X1[6][1]), p_0 = 0\$")
plot!(QQ1[:,2])
plot!(QQ1[:,3], legend=false)

QQ3 = 1 .- Q3 
P3 = plot(QQ3[:,1], xlabel="generation / 10", ylabel="\$p\$",
          title="\$m/s=$(X1[6][1]), p_0=1\$")
plot!(QQ3[:,2])
plot!(QQ3[:,3], legend=false)

M = HapDipMainlandIsland(N=N, k=k, m=0.2*s, u=u, arch=arch)
_, Q1b  = simulate(M, 220000, pinit=ones(L), drop=20000, thin=20) 

fixedpointit(M, [1.])

P1 = plot(title="\$L=$L, Ls=$Ls, N_es=$Ns, h=$h\$")
plot!(first.(X3), 1 .- getindex.(X3,2), label="\$\\mathbb{E}[p]_+\$", alpha=0.5, lw=3)
plot!(first.(X3), 1 .- getindex.(X3,3), label="\$\\mathbb{E}[p]_-\$", alpha=0.5, lw=3,
     marker=false,  xlabel="\$m/s\$", ylabel="\$p\$")
x, y, x_, y_ = findroots_ms(-s*h, -s + 2s*h, L)
plot!(P1, x, y, alpha=0.4, color=:black, label="\$\\tilde{p}_1\$")
plot!(P1, x_, y_, alpha=0.4, color=:black, ls=:dot, xlim=(0,0.4), label="\$\\tilde{p}_2\$")
plot!(first.(X2), 1 .- getindex.(X2,2), label="MRF", color=:black, alpha=0.4, marker=true, ms=1)
scatter!(first.(X1), markerstrokecolor=:black,
             1 .- map(X->mean(X[1500:end,:]), getindex.(X1,3)), color=:black,
             label="IBM(\$+\$)", legend=:bottomleft, ms=4)
scatter!(P1, [0.20, 0.25], [1 - mean(Q1b), 1- mean(Q1)], color=:orange, ms=4,
         markerstrokecolor=:orange, label="IBM(\$-\$)")

plot(P1,P2,P3,size=(900,230), layout=(1,3),titlefont=9,
     bottom_margin=5Plots.mm, left_margin=3Plots.mm)

savefig("notebooks/img/h1.svg")
savefig("/home/arthur_z/vimwiki/build/img/2023-03-29/h1.svg")

map(0.20:0.02:0.3) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
    ps = map(0:0.02:1) do p0
        1 .- fixedpointit(M, [p0])[:,:,1]
    end
    ll = maximum(length.(ps))
    P = plot()
    map(ps) do p1
        plot!([p1 ; [p1[end] for i=1:(ll-length(p1))]], title="\$m/s=$ms\$",
             color=:black, alpha=0.2)
    end
    P
end |> x->plot(x..., size=(600,300), legend=false, xlabel="iteration", ylabel="\$p\$")

savefig("notebooks/img/h1fp.svg")
savefig("/home/arthur_z/vimwiki/build/img/2023-03-29/h1fp.svg")


#
Ne2N(Ne, k) = ceil(Int, Ne/(2k) + Ne)
k  = 5
h  = 1. 
t  = 0.2
Ns = 4.  
s  = 0.02
N  = Ne2N(Ns/s, k) 
Ls = 1.5
L  = ceil(Int, Ls/s)
u  = 0.005*s
ngen = 305000
drop = 5000
thin = 50
arch = Architecture([HapDipLocus(-s*(1-t), -s*h*t, -s*t) for i=1:L])
GS(L) = GibbsSampler([UnitIntervalProposal() for i=1:L])
mmax = 1.2

X2 = tmap(0:0.1:mmax) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
    Q2 = gibbs(M, GS(L), rand(length(arch)), 10500, drop=500)
    q2 = mean(Q2)
    ss = vec(ess(Chains(Q2))[:,2]) 
    ms, q2, Q2, ss
end

X3 = tmap(0:0.01:mmax) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
    q3 = fixedpointit(M, [1.0])[end,:,:]
    q4 = fixedpointit(M, [0.0])[end,:,:]
    ms, q3[1], q4[1], q3[2], q4[2]
end

P1 = plot(title="\$L=$L, Ls=$Ls, N_es=$Ns, h=$h\$")
plot!(first.(X3), 1 .- getindex.(X3,2), label="\$\\mathbb{E}[p]_+\$", alpha=0.5, lw=3)
plot!(first.(X3), 1 .- getindex.(X3,3), label="\$\\mathbb{E}[p]_-\$", alpha=0.5, lw=3,
     marker=false,  xlabel="\$m/s\$", ylabel="\$p\$")
x, y, x_, y_ = findroots_ms(-s*(1-t) -s*h*t, -s*t + 2s*h*t, L)
plot!(P1, x, y, alpha=0.4, color=:black, label="\$\\tilde{p}_1\$")
plot!(P1, x_, y_, alpha=0.4, color=:black, ls=:dot, xlim=(0,mmax), label="\$\\tilde{p}_2\$")
plot!(first.(X2), 1 .- getindex.(X2,2), label="MRF", color=:black, marker=true, ms=2)

plot(P1, size=(350,270))
