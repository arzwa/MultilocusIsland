
using ThreadTools, Printf, Plots, PlotThemes, Serialization, MCMCChains
using Base.Iterators, Parameters, Random, StatsBase, MultilocusIsland
using ColorSchemes, Distributions, DataFrames, CSV
theme(:hokusai)

# Contrast recessive and dominant variants in the same barrier, for different
# `Ls`.
N   = 500
k   = 5
Ne  = harmonicmean(N, 2N*k)
Nes = 10  # per locus
s   = Nes/Ne
mss = 0:0.05:1
u   = s*0.005
Ls  = [10, 20, 40, 80]
τs  = [0.2, 0.8]
n   = 5000      # Gibbs samples
ker = BetaFlipProposal(0.2, 1.0, 0.1)

df1 = tmap(Iterators.product(Ls, τs, mss)) do (L, τ, ms)
    nrep = max(2, 40÷L)
    map(1:nrep) do j
        hs = [1., 0.]
        A1 = [HapDipLocus(-s*(1-τ), -s*hs[1]*τ, -s*τ) for i=1:(L÷2)]
        A2 = [HapDipLocus(-s*(1-τ), -s*hs[2]*τ, -s*τ) for i=1:(L÷2)]
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=[A1;A2])
        Q, _ = gibbs(M, ker, rand(L), n+100, drop=100)
        @info "(done) $L, $τ, $ms, $j"
        (L=L, τ=τ, h=hs, ms=ms, rep=j, 
         qs=vec(mean(Q, dims=1)), Q=cor(Q), 
         s=s, Ne=Ne, u=u, a=a, b=b)
    end
end |> x->vcat(x...) |> DataFrame

cs = ColorSchemes.Hokusai3
ps = Dict(t=>plot(grid=false) for t in τs)
xs = map(collect(groupby(df1, [:τ, :L]))) do sdf
    t, L = sdf[1,:τ], sdf[1,:L]
    c = get(cs, (L-minimum(Ls))/(maximum(Ls)-minimum(Ls)))
    xdf = combine(groupby(sdf, :ms), :qs=>ByRow(x->mean(x[1:L÷2])))
    plot!(ps[t], xdf[:,:ms], 1 .- xdf[:,2], color=c,
#          marker=true, markerstrokecolor=c, ms=2,
          label="\$L=$L\$")
    xdf = combine(groupby(sdf, :ms), :qs=>ByRow(x->mean(x[(L÷2)+1:end])))
    plot!(ps[t], xdf[:,:ms], 1 .- xdf[:,2], ls=:dot, color=c,
 #         marker=true, markerstrokecolor=c, ms=2,
          label="")
end
map(t->annotate!(ps[t], (0, 0.1, text("\$\\tau = $t\$", 7, :left))), τs)
xlabel!(ps[τs[end]], "\$m/s\$")
title!(ps[τs[1]], "\$N_es = $Nes, s=$(@sprintf "%.3f" s), \\alpha=$a, \\beta=$b\$", titlefont=8)
plot!([ps[t] for t in τs]..., layout=(3,1), size=(250,410), margin=0.2Plots.mm, 
      legend=:outertopright, legendfont=7, ylabel="\$\\mathbb{E}[p]\$")
plot(values(ps)..., size=(600,200))

# But this is not terribly interesting. More interesting would seem to contrast
# a barrier which consists mostly of haploidly selected loci (a proportion `1-τ`)
# and `τ` diploidly selected loci: how de the diploid loci behave? Does
# dominance matter at all? We could also add a third class of loci, that are
# both haploidly and diploidly selected.

# Let us consider the following: L/2 haploidly selected loci with strength `s`.
# L/2 diploidly selected loci with dominant local adaptation (`h=0`) and the
# same selection coefficient. This is the (s,h,t)-model with t∈[0,1] and h=0
N   = 1000
k   = 5
Ne  = harmonicmean(N, 2N*k)
Nes = 20  # per locus
s   = Nes/Ne
mss = 0:0.05:1
u   = s*0.005
Ls  = [2, 10, 20, 40, 80, 160]
n   = 5000      # Gibbs samples
ker = BetaFlipProposal(0.2, 1.0, 0.1)

df2 = tmap(Iterators.product(Ls, mss)) do (L, ms)
    nrep = max(2, 40÷L)
    map(1:nrep) do j
        h  = 1.0
        A1 = [HapDipLocus(-s, 0., 0.) for i=1:(L÷2)]
        A2 = [HapDipLocus(0., -h*s, -s) for i=1:(L÷2)]
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=[A1;A2])
        Q, _ = gibbs(M, ker, rand(L), n+100, drop=100)
        @info "(done) $L, $ms, $j"
        qm = vec(mean(Q, dims=1))
        qc = cor(Q)
        (L=L, h=h, ms=ms, rep=j, qs=qm, Q=qc, s=s, Ne=Ne, u=u)
    end
end |> x->vcat(x...) |> DataFrame


cs = ColorSchemes.Hokusai3
PP = plot()
xs = map(collect(groupby(df1, :L))) do sdf
    L = sdf[1,:L]
    c = get(cs, (L-minimum(Ls))/(maximum(Ls)-minimum(Ls)))
    xdf = combine(groupby(sdf, :ms), 
                  :qs=>ByRow(x->mean(x[1:L÷2])) => :q1,
                  :qs=>ByRow(x->mean(x[(L÷2)+1:end])) => :q2)
    xxdf = combine(groupby(xdf, :ms), :q1=>mean, :q2=>mean)
    plot!(PP, xxdf[:,:ms], 1 .- xxdf[:,:q1_mean], color=c, lw=2, label="\$L=$L\$")
    plot!(PP, xxdf[:,:ms], 1 .- xxdf[:,:q2_mean], color=c, lw=2, ls=:dash, alpha=0.7, label="")
end
plot(PP, size=(400,250), xlabel="\$m/s\$", ylabel="\$\\mathbb{E}[p]\$")
