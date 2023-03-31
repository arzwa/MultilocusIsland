using ThreadTools, Printf, Plots, PlotThemes, Serialization, MCMCChains
using Parameters, Random, StatsBase, MultilocusIsland
theme(:hokusai)

# Get simulation results for different dominance coefficients (`hs`) and strength
# of haploid selection (`ts`) in a parameter regime where LD matters (`Ls ~ 1`)
N = 200
k = 5
Ne = 1/(1/N + 1/(2N*k))
hs = [0.0, 0.5, 1.0]
ts = 0:0.25:1
L  = 40 
Ns = 4.  
Ls = Ns*L/Ne
mss1 = 0:0.1:1.0
mss2 = 0:0.05:1.0
mss3 = 0:0.01:1.0

# Gibbs sampler settings, for this problem these seem to be enough to get
# decent ESS
n = 5000
GS(L) = GibbsSampler([UnitIntervalProposal() for i=1:L])

# do simulations
Xs = map(hs) do h
    res = map(ts) do t
        @info (t, h)
        s = Ns/Ne
        u = s*0.01
        arch = Architecture([HapDipLocus(-s*(1-t), -s*h*t, -s*t) for i=1:L])
        # individual-based
        X = tmap(mss1) do ms
            M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
            _, Q1  = simulate(M, 11000, drop=1000, thin=10) 
            q1 = mean(Q1)
            (ms, q1)
        end
        # Gibbs
        Y = tmap(mss2) do ms
            M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
            Q2 = gibbs(M, GS(L), rand(length(arch)), n+500, drop=500)
            q2 = mean(Q2)
            ss = vec(ess(Chains(Q2))[:,2]) 
            (ms, q2, Q2, ss)
        end
        # quadrature
        Z = tmap(mss3) do ms
            M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
            q3 = fixedpointit(M, [0.0])[end,1,1]
            (ms, q3)
        end
        (t, X, Y, Z)
    end
    (h, res)
end

# Save the simulations
serialize("data/dominance-haploidsel.jls", Xs)

# Make a plots
map(enumerate(Xs)) do (j,(h, Y))
    P = plot()
    map(enumerate(Y)) do (i,(t, X, Y, Z))
        c = i
        Z = Ys[j][2][i][2]
        plot!(first.(Z), 1 .- last.(Z), 
              color=c, title="\$h=$h\$",#, L=$L, N_es=$Ns\$", 
              label="",lw=2.5, alpha=0.3)
        gs = map(getindex.(Y, 3)) do xg
        end
        plot!(first.(Y), 1 .- getindex.(Y,2), 
              ls=:dot, color=c, label="")
        scatter!(first.(X), 1 .- getindex.(X, 2),
                 label="\$\\tau=$(@sprintf "%.2f" t)\$",
                 color=c, markerstrokecolor=c,
                 xlabel="\$m/s\$", ylabel="\$\\mathbb{E}[p]\$")
    end 
    P
end |> x->plot(x..., grid=false, ms=3, ylim=(0,1), xlim=(0,1.05), 
               margin=3Plots.mm, layout=(1,3), legend=:topright,
              size=(700,200), legendfont=6)

#savefig("img/domtauL40Nes4.pdf")
savefig("/home/arthur_z/vimwiki/build/img/2023-03-30/domtau.svg")


# fixed point iteration starting from no differentiation
Ys = map(hs) do h
    res = map(ts) do t
        @info (t, h)
        s = Ns/Ne
        u = s*0.01
        arch = Architecture([HapDipLocus(-s*(1-t), -s*h*t, -s*t) for i=1:L])
        Z = tmap(mss3) do ms
            M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
            q3 = fixedpointit(M, [0.0])[end,1,1]
            (ms, q3)
        end
        (t, Z)
    end
    (h, res)
end

Xs = deserialize("data/dominance-haploidsel.jls")

map(enumerate(Xs)) do (j,(h, Y))
    P = plot()
    map(enumerate(Y)) do (i,(t, X, Y, Z))
        c = i
        plot!(first.(Z), 1 .- last.(Z), 
              color=c, title="\$h=$h\$",#, L=$L, N_es=$Ns\$", 
              label="",lw=2.5, alpha=0.3)
        ZZ = Ys[j][2][i][2]
        plot!(first.(ZZ), 1 .- last.(ZZ), 
              color=c, title="\$h=$h\$",#, L=$L, N_es=$Ns\$", 
              label="",lw=2.5, alpha=0.3)
        plot!(first.(Y), 1 .- getindex.(Y,2), 
              ls=:dot, color=c, label="")
        scatter!(first.(X), 1 .- getindex.(X, 2),
                 label="\$\\tau=$(@sprintf "%.2f" t)\$",
                 color=c, markerstrokecolor=c,
                 xlabel="\$m/s\$", ylabel="\$\\mathbb{E}[p]\$")
    end 
    P
end |> x->plot(x..., grid=false, ms=3, ylim=(0,1), xlim=(0,1.05), 
               margin=3Plots.mm, layout=(1,3), legend=:topright,
              size=(700,200), legendfont=6)


k = 5
Ne = 1/(1/N + 1/(2N*k))
h  = 0.9 
t  = 1.
Ns = 4.  
Ls = 0.88
mss1 = 0:0.1:1.0
mss2 = 0:0.05:1.0
mss3 = 0:0.01:1.0
Ne2N(Ne, k) = ceil(Int, Ne/(2k) + Ne)

n = 20000
GS(L) = GibbsSampler([UnitIntervalProposal() for i=1:L])

XX = map([20,40,60,80]) do L
    s = Ls/L
    N = Ne2N(Ns/s, k)
    arch = Architecture([HapDipLocus(-s*(1-t), -s*h*t, -s*t) for i=1:L])
    @show (L, s, N)
    u = s*0.01
    X = tmap(mss1) do ms
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
        _, Q1  = simulate(M, 21000, drop=1000, thin=10) 
        q1 = mean(Q1)
        (ms, q1)
    end
    Y = tmap(mss2) do ms
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
        Q2, _ = gibbs(M, GS(L), rand(length(arch)), n+100, drop=100)
        q2 = mean(Q2)
        ss = vec(ess(Chains(Q2))[:,2]) 
        (ms, q2, Q2, ss)
    end
    # quadrature
    Z = tmap(mss3) do ms
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
        q3 = fixedpointit(M, [1.0])[end]
        (ms, q3)
    end
    L, X, Y, Z
end

plot()
map(enumerate(XX)) do (i,(L, X, Y, Z))
    plot!(first.(Y), getindex.(Y, 2), color=i, ls=:dot, label="")
    plot!(first.(Z), last.(Z), color=i, alpha=0.3, lw=2, label="\$L=$L\$")
    scatter!(first.(X), last.(X), color=i, label="", markerstrokecolor=i, ms=4)
end
plot!(xlabel="\$m/s\$", ylabel="\$q\$", size=(300,230), legend=:bottomright)


savefig("/home/arthur_z/vimwiki/build/img/2023-03-27/h1fp.pdf")
