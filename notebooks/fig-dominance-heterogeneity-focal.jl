using MultilocusIsland, StatsBase, ThreadTools, Printf, QuadGK, StatsPlots, Serialization
using Plots, PlotThemes, Distributions, Parameters; theme(:hokusai)

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

map(enumerate(results)) do (i,xh)
    P = plot(legend=i == 3 ? :right : false)
    map(enumerate(xh)) do (j,xs)
        h,s,ys1,ys2 = xs
        map(zip(ys1,mss1)) do (xs,ms)
            scatter!(repeat([ms], length(xs)), xs, alpha=0.2, ms=2, color=j,
                     label="")
        end
        plot!(mss1, map(mean, ys1), color=j, label="", title="\$h=$h\$")
        plot!(mss2, ys2, color=j, alpha=0.5, lw=4, label="\$s=$(@sprintf "%.3f" s)\$")
    end
    P
end |> x->plot(x..., layout=(1,3), size=(750,200), xlabel="\$m/\\bar{s}\$",
               ylabel="\$\\mathbb{E}[p]\$", ylim=(0,1), xlim=(0,1.25), margin=3Plots.mm)


function randbarrier(df, L, mx=0.3)
    ss = rand(truncated(df.sd, 0.0, mx), L)
    hs = rand(df.hd, L)
    ss = ss .* mean(df.sd)/(mean(ss))
    hs = hs .* mean(df.hd)/(mean(hs))
    [DipLocus(-s*h, -s) for (s,h) in zip(ss,hs)]
end

# Should compare homogeneous vs heterogeneous barrier effect on a focal
# locus... Didn't do this for some reason?
κs = [4, 1, 1/4]
ss = [s̄/2, s̄, s̄*4]
hs = [0.0, 0.5, 1.0]
mss1 = Dict(ss[1] => 0.05:0.10:1.05, ss[2]=>0.05:0.1:1.05, ss[3]=> 0.1:0.2:2.1)
nrep = 100

results = map(hs) do h
    map(ss) do s
        map(κs) do κ
            sd = Gamma(κ, s̄/κ)
            df = IndependentDFE(sd, Beta(1, 1))
            @info κ, s, h
            ys1 = map(mss1[s]) do ms
                @info ms
                # random architectures
                xs_ = tmap(1:nrep) do _
                    A  = Architecture([randbarrier(df, L-1) ; DipLocus(-s*h,-s)])
                    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), ms*s̄)
                    try 
                        P,_= fixedpointit(M, ones(L));
                        pm = P[end,:,1]
                    catch
                        pm = [missing for i=1:L]
                    end
                end
                xs = filter(x->!all(ismissing.(x)), xs_)
                hcat(xs...)
            end
            (κ, h, s, ys1)
        end |> x->vcat(x...)
    end
end

map(enumerate(results)) do (i,hx)
    y = map(enumerate(hx)) do (j,sx)
        P = plot(legend=j==1 ? :topright : false)
        map(enumerate(sx)) do (kk,ks)
            κ,h,s,ys1 = ks
            k = [1,3,5][kk]
            zs = map(enumerate(zip(ys1,mss1[s]))) do (m,(xs,ms))
                ps = xs[end,:]
                violin!([m], ps, color=k, label="", width=0.1, alpha=0.5, linealpha=0.0, linecolor=k)
            end
            scatter!(map(y->mean(y[end,:]), ys1), color=k, label="\$\\kappa=$κ\$",
                  title="\$h=$h, s=$s\$", ylim=(0,1))#, xlim=(0,maximum(mss1[s])+0.05))
            plot!(map(mean, ys1), color=k, label="", lw=2, alpha=1, marker=true, ms=2,
                  title="\$h=$h, s=$s\$", ylim=(0,1))#, xlim=(0,maximum(mss1[s])+0.05))
            xticks!(1:2:length(mss1[s]), string.(mss1[s][1:2:end]))
        end
        # single locus prediction
        ms_ = range(extrema(mss1[s])..., 200)
        xs_ = range(1, length(mss1[s]) , 200)
        c = length(ms_)/length(mss1[s])
        _,h,s,_ = sx[1]
        zs  = map(ms_) do ms
            A = Architecture([DipLocus(-s*h, -s)])
            M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), ms*s̄)
            Q,_= fixedpointit(M, ones(1))
            Q[end,1,1]
        end
        plot!(xs_, zs, color=:black, alpha=0.2, lw=2, label="")
        P
    end |> x->vcat(x...)
end |> x->plot(vcat(x...)..., layout=(3,3), size=(700,600),
               ylabel="\$\\mathbb{E}[p]\$", xlabel="\$m/\\bar{s}\$")

results = deserialize("data/hetfocal3.jls")


# Heterogeneity on the x-axis
# ===========================
ks = collect(3:-0.5:-3)
κs = 2 .^ ks
nrep = 100
ms = 0.2:0.2:1.0
s = s̄
hs = [0.0, 0.5, 1.0]

out = map(hs) do h
    map(ms) do m
        map(κs) do κ
            @info (h, m, κ)
            sd = Gamma(κ, s̄/κ)
            df = IndependentDFE(sd, Beta(1, 1))
            xs_ = tmap(1:nrep) do _
                A  = Architecture([randbarrier(df, L-1) ; DipLocus(-s*h,-s)])
                M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄)
                try 
                    P,_= fixedpointit(M, ones(L));
                    pm = P[end,:,1]
                catch
                    pm = [missing for i=1:L]
                end
            end
            xs = filter(x->!all(ismissing.(x)), xs_)
            hcat(xs...)
        end
    end
end

map(zip(hs,out)) do (h, Y)
    plot(title="\$h=$h, s=\\bar{s}=$s̄, L\\bar{s}=$Ls, N_e\\bar{s}=$Ns\$")
    map(enumerate(zip(ms, Y))) do (i,(m, X))
        map(enumerate(zip(κs, X))) do (j,(k,xs))
            violin!([j], xs[end,:], color=i, alpha=0.5, linecolor=i, linealpha=0.0,
                    label="")
        end
        plot!(map(x->mean(x[end,:]), X), marker=true, color=i, ms=4, lw=2,
              label="\$m/\\bar{s} = $m\$")
        # single locus prediction
        A = Architecture([DipLocus(-s*h, -s)])
        M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄)
        P,_= fixedpointit(M, ones(1))
        hline!([P[end,1,1]], color=i, label="", alpha=0.4)
    end
    plot!(xticks=(1:2:length(κs), [@sprintf("\$2^{%d}\$", k) for k in ks[1:2:end]]), 
          xlabel="\$\\kappa\$", ylabel="\$\\mathbb{E}[p]\$", ylim=(0,1), size=(500,300))
end |> x->plot(x..., layout=(3,1), size=(350,700), margin=3Plots.mm)

out = deserialize("data/hetfocalb.jls")
