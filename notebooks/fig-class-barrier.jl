using MultilocusIsland
using Plots, PlotThemes, Printf, ColorSchemes; theme(:hokusai)
cs = ColorSchemes.viridis

s = 0.02
l = 15
h = [zeros(l); fill(0.5,l); ones(l)]
L = length(h)
τ = 1.
N = 5000
k = 5
u = s*0.005
m = 0.01

A1 = Architecture([HapDipLocus(-s*(1-τ), -s*h[i]*τ, -s*τ) for i=1:L], fill(0.5, L)) 
A2 = Architecture([HapDipLocus(-s*(1-τ), 0., -s*τ)      for i=1:L], fill(0.5, L)) 
A3 = Architecture([HapDipLocus(-s*(1-τ), -s*τ*0.5, -s*τ)  for i=1:L], fill(0.5, L)) 
A4 = Architecture([HapDipLocus(-s*(1-τ), -s*τ, -s*τ)        for i=1:L], fill(0.5, L)) 
    
mss = 0:0.005:1
ps = map(mss) do ms
    xs = map(zip([A1,A2,A3,A4], [3,1,1,1])) do (A,K)
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, arch=A, u=u)
        P,_ = fixedpointit(M, ones(K))
        P[end,:,1]
    end |> x->vcat(x...)
end |> x->hcat(x...) |> permutedims 

plot( mss, ps[:,1:3], label=["dominant" "codominant" "recessive"])
plot!(mss, ps[:,4:6], color=[1 2 3], ls=:dot, label="")
plot!(mss, mean(ps[:,1:3], dims=2), color=:black, lw=2, label="average")
plot!(xlabel="\$m/s\$", ylabel="\$p\$")

# large effect and a bunch of small effect loci
ss = [fill(0.01, 35); fill(0.05, 5)] 
A1 = Architecture([HapDipLocus(-s*(1-τ), -s*τ, -s*τ) for s in ss]) 

mss = 0:0.005:1
ps = map(mss) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, arch=A1, u=u)
    fixedpointit(M, ones(2))[end,:,1]
end |> x->hcat(x...) |> permutedims 
plot(ps)


# a single recessive, a single dominant, and a bunch of additives
Ls = 1.
L  = 50
s  = Ls/L
Ns = 10.
Ne = Ns/s
k  = 5
N  = _Ne2N(Ne, k) 

P1s = map([1, 0]) do h
    A  = Architecture(DipLocus(-0.5s, -s), L) 
    push!(A, DipLocus(-s*h, -s), 0.5)
    mss = 0:0.005:0.65
    ps = map(mss) do ms
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, arch=A, u=s*0.005)
        P,_ = fixedpointit(M, ones(2))
        P[end,:,1]
    end |> x->hcat(x...) |> permutedims 
    lab = h == 0. ? "dominant" : (h == 1 ? "recessive" : "\$h = $h\$")
    P1 = plot(mss, ps, color=[2 4], 
              label=["codominant" lab], legend=:topright,
              xlabel="\$m/s\$", ylabel="\$\\mathbb{E}[p]\$", lw=2, 
              title="\$N_es = $Ns, Ls=$Ls, L=$L\$")
    ps2 = map(mss) do ms
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, arch=A[L+1:L+1], u=s*0.005)
        P,_= fixedpointit(M, ones(1));
        P[end,1,1]
    end
    plot!(P1, mss, ps2, color=4, label="single locus", ls=:dash)
end

P2s = map([1, 0]) do h
    lab = h == 0. ? "dominant" : (h == 1 ? "recessive" : "\$h = $h\$")
    A  = Architecture(DipLocus(-0.5s, -s), L) 
    push!(A, DipLocus(-s*h, -s), 0.5)
    sm = 5
    m = s/sm
    M = HapDipMainlandIsland(N=N, k=k, m=m, arch=A, u=s*0.005)
    P,_= fixedpointit(M, ones(2));
    p  = P[end,:,1]
    pq = P[end,:,2]
    ys = expectedsfs(M, p, pq, step=0.005, f=log10)
    M = HapDipMainlandIsland(N=N, k=k, m=m, arch=A[L+1:L+1], u=s*0.005)
    P,_= fixedpointit(M, ones(1));
    p  = P[end,:,1]
    pq = P[end,:,2]
    ys2 = expectedsfs(M, p, pq, step=0.005, f=log10)
    P2 = plot(ys[2], color=:black, lw=2, label="multilocus", title="\$s = $(sm)m\$, $lab")
    plot!(ys2, color=:black, alpha=0.4, lw=2, label="single locus", xlabel="\$p\$",
          ylabel="\$\\log_{10}\\phi(p)\$",
          size=(300,250), legend=:bottomright)
end

plot(P1s[1], P2s[1], size=(500,200), margin=3Plots.mm)
plot(P1s[2], P2s[2], size=(500,200), margin=3Plots.mm)


