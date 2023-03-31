using ThreadTools, Printf, Plots, PlotThemes, Serialization, MCMCChains
using Base.Iterators, Parameters, Random, StatsBase, MultilocusIsland
using ColorSchemes, Distributions, DataFrames, CSV
theme(:hokusai)

# We want to sample appropriate dominance coefficients for a given τ. There are
# obviously multiple ways to do this.

# This function assumes a single shared selection coefficient, and a branching
# process approximation for the fixation probability. Under these assumptions
# the distribution of the fixed `h` does not depend on `s` nor `Ne`. It does
# give a pronounced Haldane's sieve effect (much more so than when we use a
# diffusion-based fixation probability for not very large populations).
# Note that the distribution has a closed form, but it is more straightforward
# to sample from a sufficiently fine discretization.
function sample_h1(hdist, τ, n)
    ps = map(h->(1-τ*(1-h))*pdf(hdist, h), 0:0.01:1)
    ps ./= sum(ps)
    1 .- rand(DiscreteNonParametric(0:0.01:1, ps), n)
end

# **Note** that the dominance coefficients `h` for beneficials fixed on the
# island in allopatry correspond to dominance coefficients `1-h` of invading
# (ancestral) alleles from the mainland in the mainland-island model.

# The expected value of the dominance coefficient
function expectedh(hdist::Beta, t)
    a = hdist.α
    b = hdist.β
    Eh = a * (a + 1 + b*(1-t)) / ((a + b*(1-t)) * (a + b + 1))
    return 1 - Eh
end

# Let us start by contrasting a haploid and diploid life cycle (`τ=0` vs.
# 'τ=1`). We do the computation with variable `h` and a shared `h` set to the
# expected value.
N   = 500
k   = 5
Ne  = harmonicmean(N, 2N*k)
Nes = 5  # per locus
s   = Nes/Ne
mss = 0:0.05:1
u   = s*0.005
Ls  = [1, 10, 20, 40, 80, 160]
τs  = [0., 0.5, 1.0]
n   = 5000      # Gibbs samples
a,b = 1, 1
hd  = Beta(a,b) # distribution of `h` for new mutations
GS(L) = GibbsSampler([UnitIntervalProposal() for i=1:L])

df1 = tmap(Iterators.product(Ls, τs, mss)) do (L, τ, ms)
    nrep = max(2, 40÷L)
    map(1:nrep) do j
        h = sample_h1(hd, τ, L)
        A = Architecture([HapDipLocus(-s*(1-τ), -s*h[i]*τ, -s*τ) for i=1:L])
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A)
        Q = gibbs(M, GS(L), rand(L), n+100, drop=100)
        @info "(done) $L, $τ, $ms, $j"
        (L=L, τ=τ, h=h, ms=ms, rep=j, 
         q=mean(Q), qs=vec(mean(Q, dims=1)), Q=cor(Q), 
         s=s, Ne=Ne, u=u, a=a, b=b)
    end
end |> x->vcat(x...) |> DataFrame

df2 = tmap(Iterators.product(Ls, τs, mss)) do (L, τ, ms)
    nrep = max(2, 40÷L)
    map(1:nrep) do j
        h = fill(expectedh(hd, τ), L)
        A = [HapDipLocus(-s*(1-τ), -s*h[i]*τ, -s*τ) for i=1:L]
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A)
        Q, _ = gibbs(M, ker, rand(L), n+100, drop=100)
        @info "(done) $L, $τ, $ms, $j"
        (L=L, τ=τ, h=h, ms=ms, rep=j, 
         q=mean(Q), qs=vec(mean(Q, dims=1)), Q=cor(Q), 
         s=s, Ne=Ne, u=u, a=a, b=b)
    end
end |> x->vcat(x...) |> DataFrame

cs = ColorSchemes.Hokusai3
ps = Dict(t=>plot(grid=false) for t in τs)
xs = map(collect(groupby(df1, [:τ, :L]))) do sdf
    t, L = sdf[1,:τ], sdf[1,:L]
    xdf = combine(groupby(sdf, :ms), :q=>mean)
    c = get(cs, L/maximum(Ls))
    plot!(ps[t], xdf[:,:ms], 1 .- xdf[:,2], 
          marker=true, markerstrokecolor=c, color=c, ms=2,
          label="\$L=$L\$")
end
xs = map(collect(groupby(df2, [:τ, :L]))) do sdf
    t, L = sdf[1,:τ], sdf[1,:L]
    xdf = combine(groupby(sdf, :ms), :q=>mean)
    c = get(cs, L/maximum(Ls))
    plot!(ps[t], xdf[:,:ms], 1 .- xdf[:,2], alpha=0.2, label="")
end
map(t->annotate!(ps[t], (0, 0.1, text("\$\\tau = $t\$", 7, :left))), τs)
xlabel!(ps[τs[end]], "\$m/s\$")
title!(ps[τs[1]], "\$N_es = $Nes, s=$(@sprintf "%.3f" s), \\alpha=$a, \\beta=$b\$", titlefont=8)
plot!([ps[t] for t in τs]..., layout=(3,1), size=(250,410), margin=0.2Plots.mm, 
      legend=:outertopright, legendfont=7, ylabel="\$\\mathbb{E}[p]\$")

savefig("notebooks/img/dfe-divergence-1.pdf")
savefig("/home/arthur_z/vimwiki/notes/img/dfe-divergence-1b.svg")

CSV.write("data/dfe-divergence-1.csv", df[:,Not(:Q)])


# Small differences in τ, do they matter?
N   = 500
k   = 5
Ne  = harmonicmean(N, 2N*k)
Nes = 5  # per locus
s   = Nes/Ne
mss = 0:0.025:1
u   = s*0.005
Ls  = [20, 40, 80, 160]
τs  = [0.1, 0.2]
n   = 10000      # Gibbs samples
rep = 1:3        # number of simulation replicates
a,b = 1, 3
hd  = Beta(a,b)  # distribution of `h` for new mutations
ker = BetaProposal(0.2)

Xs = tmap(Iterators.product(Ls, τs, mss, rep)) do (L, τ, ms, j)
    @info L, τ, ms, j
    h = sample_h1(hd, τ, L)
    A = [HapDipLocus(-s*(1-τ), -s*h[i]*τ, -s*τ) for i=1:L]
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A)
    Q, _ = gibbs(M, ker, rand(L), n+100)
    Q = Q[101:end,:]
    @info "(done)"
    (L=L, τ=τ, h=h, ms=ms, rep=j, q=mean(Q), Q=Q)
end

df = DataFrame(Xs)
plot(grid=false)
cs = ColorSchemes.Hokusai3
pp = plot()
xs = map(collect(groupby(df, [:τ, :L]))) do sdf
    t, L = sdf[1,:τ], sdf[1,:L]
    xdf = combine(groupby(sdf, :ms), :q=>mean)
    plot!(pp, xdf[:,:ms], 1 .- xdf[:,2], 
         label=t==0.1 ? "\$L=$L\$" : "",  
         ls = t == 0.1 ? :solid : :dot,
         color=get(cs, (L-minimum(Ls))/(maximum(Ls) - minimum(Ls))))
end
title!(pp, "\$N_es = $Nes, s=$(@sprintf "%.3f" s), \\alpha=$a, \\beta=$b\$",
       titlefont=8)
plot!(pp, size=(250,200), legend=:right, legendfont=7, 
      ylabel="\$\\mathbb{E}[p]\$", xlabel="\$m/s\$")

savefig("notebooks/img/dfe-divergence-3.pdf")
savefig("/home/arthur_z/vimwiki/build/img/dfe-divergence-3.svg")
