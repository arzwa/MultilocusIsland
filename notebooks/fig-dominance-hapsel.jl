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
n = 10000
kernel = BetaProposal(0.3)

# do simulations
Xs = map(hs) do h
    res = map(ts) do t
        @info (t, h)
        s = Ns/Ne
        u = s*0.01
        arch = [HapDipLocus(-s*(1-t), -s*h*t, -s*t) for i=1:L]
        # individual-based
        X = tmap(mss1) do ms
            M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
            _, Q1  = simulate(M, 11000) 
            q1 = mean(Q1[1000:10:end,:])
            (ms, q1)
        end
        # Gibbs
        Y = tmap(mss2) do ms
            M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
            Q2, _ = gibbs(M, kernel, rand(length(arch)), n+100)
            q2 = mean(Q2[101:end,:])
            ss = vec(ess(Chains(Q2[101:end,:]))[:,2]) 
            (ms, q2, Q2, ss)
        end
        # quadrature
        Z = tmap(mss3) do ms
            M = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=u, arch=arch)
            q3 = UnlinkedHapDip.expectedq(M, [0.01, 0.5, 0.99])[1]
            (ms, q3)
        end
        (t, X, Y, Z)
    end
    (h, res)
end

# Save the simulations
serialize("data/dominance-haploidsel.jls")

# Make a plots
map(enumerate(Xs)) do (i,(h, Y))
    P = plot()
    map(enumerate(Y)) do (i,(t, X, Y, Z))
        c = i
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

