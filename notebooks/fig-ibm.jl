# Comparisons against individual based simulations
# - assess dominance, Ls, Ns effects
# - consider arbitrary variable barrier

using ThreadTools, Printf, Plots, PlotThemes, Serialization, MCMCChains
using Base.Iterators, Parameters, Random, StatsBase, MultilocusIsland
using ColorSchemes, Distributions, DataFrames, CSV
theme(:hokusai)

function dosims(xs, mss1, mss2, mss3, kernel, ns, ng; init=[0.01, 0.5, 0.99])
    df = tmap(xs) do x
        @unpack L, Ne, N, k, s, h, τ, u = x
        A = [HapDipLocus(-s*(1-τ), -s*h*τ, -s*τ) for i=1:L]
        q1 = map(mss1) do ms
            M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A)
            _, Q1  = simulate(M, ng+1000, drop=1000, thin=10) 
            q = 1 - mean(Q1)
            (ms, q)
        end
        q2 = map(mss2) do ms
            M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A)
            Q, _ = gibbs(M, kernel, rand(L), ns+1000, drop=1000)
            q = 1 - mean(Q)
            (ms, q)
        end
#        q3 = map(mss3) do ms
#            M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A)
#            q = expectedq(M, init)[1]
#            (ms, 1 - q)
#        end
        @info "done $x"
        #(; q1=q1, q2=q2, q3=q3, x...)
        (; q1=q1, q2=q2, x...)
    end |> DataFrame
end

# 1. Check what happens as `N` gets small, keeping Ls appreciable (here 2).
# We expect the diffusion approximation to break down (N gets small).
L  = 40
s  = 2.0./L
u  = s * 0.05
Nn = [100, 50, 25, 10] 
hs = [0.0, 1.0]
k  = 5
τ  = 0.5
xs = map(Iterators.product(Nn, hs)) do (N, h)
    Ne = harmonicmean(N, 2N*k)
    (L=L, Ne=Ne, Ns=N*s, Ls=L*s, N=N, k=k, s=s, h=h, τ=τ, u=u, n=5)
end |> vec

mss1 = 0:0.10:1   # individual-based
mss2 = 0:0.05:1   # gibbs sampler
mss3 = 0:0.02:1   # numerical expectation
kernel = BetaProposal(0.2)
df = dosims(xs, mss1, mss2, mss3, kernel, 5000, 5000, init=[0.01, 0.1, 0.5, 0.9, 0.99])

map(collect(groupby(df, :N))) do sdf
    N = sdf[1,:N]
    P = plot(title="\$N=$N, L=$L, Ls=$(s*L)\$", titlefont=8, 
             legend= N==10 ? :topright : false)
    for (i,row) in enumerate(eachrow(sdf))
        c = row[:h] == 0. ? 5 : 2
        scatter!(row[:q1], ms=3, markerstrokecolor=c, color=c, label="\$h=$(row[:h])\$")
        plot!(row[:q2], color=c, label="")
#        plot!(row[:q3], color=c, ls=:dash, label="")
    end
    P
end |> x->plot(x..., size=(400,350), ylim=(0,1.), 
               ylabel="\$\\mathbb{E}[p]\$", xlabel="\$m/s\$")

savefig("/home/arthur_z/vimwiki/notes/img/driftbreakdown.pdf")

# 2. Check adequacy for increasing Ls for different dominance (in diploids).
# We include over and underdominance.
s   = 0.05
LL  = [10,20,40,80]
hs  = [-0.5, 0., 0.5, 1.0, 1.5]
us  = 0.005
k   = 5
τ   = 1.0
N   = 200
Ne  = harmonicmean(N, 2N*k)
nrep= 80

xs = map(Iterators.product(LL, hs)) do (L, h)
    n = max(1, nrep÷L)
    u = us*s
    [(L=L, Ne=Ne, Ns=Ne*s, Ls=L*s, N=N, k=k, s=s, h=h, τ=τ, u=u) for i=1:n]
end |> x->mapreduce(y->vcat(y...), vcat, x)

mss1 = 0:0.10:1   # individual-based
mss2 = 0:0.05:1   # gibbs sampler
mss3 = 0:0.50:1   # numerical expectation
kernel = BetaFlipProposal(0.5, 1.0, 0.1)
ns = 5000  # gibbs samples
ng = 5000  # IBM generations
df = dosims(xs, mss1, mss2, mss3, kernel, ns, ng, init=[0.1, 0.5, 0.9])

df4 = dosims(xs, mss1, mss2, mss3, kernel, ns, ng, init=[0.1, 0.5, 0.9])

df = vcat(df, df4)

function avgtups(xs::Vector)
    ls = unique(length.(xs))
    @assert length(ls) == 1
    xx = xs[1]
    x1 = first.(xx)
    zs = map(2:length(xx[1])) do i
        ys = map(x->getindex.(x, i), xs)
        mean(ys)
    end
    collect(zip(x1, zs...))
end

Ps = map(collect(groupby(df, :L))) do sdf
    L = sdf[1,:L]
    P = plot(title="\$L=$L, Ls=$(s*L)\$", titlefont=8, legend=L==10 ? :topright : false)
    map(collect(groupby(sdf, :h))) do ssdf
        h = ssdf[1,:h]
        c = Dict(h=>i for (i,h) in enumerate(hs))[h]
        q1 = avgtups(ssdf[:,:q1])
        q2 = avgtups(ssdf[:,:q2])
        scatter!(q1, ms=3, markerstrokecolor=c, color=c, label="\$h=$h\$")
        plot!(q2, color=c, label="")
    end
    P
end 
plot(Ps..., size=(500,380), ylim=(0,1.), ylabel="\$\\mathbb{E}[p]\$",
     xlabel="\$m/s\$")

savefig("/home/arthur_z/vimwiki/notes/img/dominance-Ls.pdf")

# Hard cases
function avgcor(Q)
    C = cor(Q)
    c = 0.
    n = size(C,1)
    for i=2:n
        for j=1:i-1
            c += C[i,j]
        end
    end
    2c/(n*(n-1))  
end

xx = map(hs) do h
    @info h
    L = 20
    s = 0.05
    u = s*0.005
    k = 5
    N = 200
    τ = 1.
    A = Architecture([HapDipLocus(-s*(1-τ), -s*h*τ, -s*τ) for i=1:L], fill(0.5, L))
    mss = 0:0.05:1
    qs = map(mss) do ms
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=A)
        #K = BetaFlipProposal(0.5, 1.0, 0.1)
        G = GibbsSampler([UnitIntervalProposal() for i=1:L])
        Q, _ = gibbs(M, G, rand(L), 11000, drop=1000)
        mean(Q), avgcor(Q), Q
    end
    q = 1 .- first.(qs)
    c = getindex.(qs, 2)
    q, c
end

Ps = map(1:5) do i
    q, c = xx[i]
    plot(0:0.05:1, q, legend=:topright,
         color=:black, xlabel="\$m/s\$", title="\$h=$(hs[i])\$",
         label="\$\\mathbb{E}[p]\$", ylim=(0,1))
    plot!(0:0.05:1, c, color=:salmon, label="\$\\rho\$")
end 
plot(Ps..., size=(500,270))

savefig("/home/arthur_z/vimwiki/notes/img/dominance-cor.pdf")

plot!(Ps[2], mss, q, color=:red)
plot(Ps[2])
