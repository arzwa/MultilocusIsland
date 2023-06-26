using MultilocusIsland, StatsBase, ColorSchemes, ThreadTools, Printf, QuadGK
using Plots, PlotThemes, Distributions, Parameters; theme(:hokusai)
using LogExpFunctions
cs = ColorSchemes.viridis

map([2/3, 1, 3/2]) do κ
    s̄ = 0.01
    λ = κ/s̄
    dfe1 = IndependentDFE(Gamma(κ, 1/λ), Beta(1, 1))
    dfe2 = Logisticsbyh(Gamma(κ, 1/λ), (s̄/2, 0.5), (0.1, 0.99), 1.)
    dfe3 = CKGamma(κ, λ, 1/3)
    ss = 0.0001:0.001:0.1
    hs = 0.0:0.01:1
    P1 = contourf(ss, hs, (s,h)->logpdf(dfe1, s, h), levels=40, lines=false, title=κ == 2/3 ? "(A)" : "")
    P2 = contourf(ss, hs, (s,h)->logpdf(dfe2, s, h), levels=40, lines=false, title=κ == 2/3 ? "(B)" : "")
    P3 = contourf(ss, hs, (s,h)->logpdf(dfe3, s, h), levels=40, lines=false, title=κ == 2/3 ? "(C)" : "")
    plot!(P3, x->1-exp(-dfe3.K*x), color=:black, lw=1, legend=false)
    P  = plot(P1, P2, P3, size=(700,200), layout=(1,3), colorbar=false, xlim=(0,0.1001), tickfont=7,
              ylim=(0,1), xlabel="\$s\$", ylabel="\$h\$")#, bottom_margin=3Plots.mm, right_margin=0Plots.mm)
end |> x->plot(x..., layout=(3,1), size=(600,500), margin=0Plots.mm)

map([2/3, 1., 3/2]) do κ
    s̄ = 0.01
    λ = κ/s̄
    dfe1 = IndependentDFE(Gamma(κ, 1/λ), Beta(1, 1))
    dfe2 = Logisticsbyh(Gamma(κ, 1/λ), (s̄/2, 0.5), (0.1, 0.99), 1.)
    dfe3 = CKGamma(κ, λ, 1/3)
    P1 = plot(dfe1.sd, color=:black, fill=true, fillalpha=0.2, legend=false, xlabel="\$s\$",
              title="\$\\bar{s} = $s̄, \\kappa = $(@sprintf "%.2f" κ)\$", xlim=(-0.003,0.103))
    P2 = plot()
    map(zip([dfe1, dfe2, dfe3], ["independent", "logistic", "CK94"])) do (df, lab)
        plot!(P2, 0:0.01:1, h->MultilocusIsland.pdfh(df, h), fill=true,
              fillalpha=0.2, label=lab, legend=:topleft, xlabel="\$h\$")
    end
    plot(P1, P2, layout=(2,1))
end |> x->plot(x..., layout=(1,3), size=(600,300))

plot( ss, s->logistic(dfe2.a + dfe2.b*log(s)))
plot!(ss, s->1-exp(-s*dfe3.K)/2)

function randsh(dfe)
    l = MultilocusIsland.randlocus(dfe)
    -l.s11, l.s01/l.s11
end

S0 = scatter(map(_->randsh(dfe1), 1:10000), color=:gray, alpha=0.3)
S1 = scatter(map(_->randsh(dfe2), 1:10000), color=:gray, alpha=0.3)
S2 = scatter(map(_->randsh(dfe3), 1:10000), color=:gray, alpha=0.3)
S  = plot(S0, S1, S2, size=(700,200), layout=(1,3), colorbar=false, xlim=(0,0.151),
          ylim=(0,1), xlabel="\$s\$", ylabel="\$h\$", legend=false, ms=1,
          right_margin=5Plots.mm)# margin=3Plots.mm, ms=1, legend=false)
plot(P, S, layout=(2,1), size=(700,400))

plot(0:0.005:1, map(h->MultilocusIsland.marginalh(dfe2, h), 0:0.005:1))

# DFE at migration-selection balance
#
# A way that might be better to characterize the DFE at migration-selection
# equilibrium: consider one samples an individual haploid genome at equilibrium
# on the island, given that a locus is differentiated, what is the probability
# distribution over `s` and `h`? this is roughly what we did above, but it can
# be made more precise by taking this operational definition serious.

# computational intensive, but more proper, way
function msedfe(dfe, ss, hs, L, N, k, u, m; n=100)
    Z = 0.
    Y = map(ss) do s
        map(hs) do h
            xs = tmap(1:n) do _
                B = [randlocus(dfe) for i=1:L-1]
                A = Architecture([DipLocus(-s*h, -s) ; B])
                M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m, ones(L))
                P,_= fixedpointit(M, ones(L));
                pback = P[end,2:end,1]
                pfore = P[end,1,1]
                Z += sum(pback)
                pfore * pdf(dfe, s, h)
            end
            sum(xs)/n
        end 
    end 
    ED = Z / (length(ss) * length(hs) * n * L)
    ED, hcat(Y...) ./ ED
end

function singlelocuseq(M)
    map(1:length(M.D.A)) do i
        DD = reconstruct(M.D, A=M.D.A[i:i])
        MM = MainlandIslandModel(DD, M.m, ones(1))
        P,_= fixedpointit(MM, ones(1))
        P[end,:,1]
    end |> x->first.(x)
end

# amounts more or les sto the same I suppose.
function msedfe1(dfe, L, N, k, u, m; n=100) 
    sss = Float64[]
    hss = Float64[]
    wss = Float64[]
    dss = Float64[]
    for i=1:n
        A  = Architecture([randlocus(dfe) for i=1:L])
        sa = [l.s01 for l in A.loci]
        sb = [l.s11 - 2l.s01 for l in A.loci]
        ss = [l.s11 for l in A.loci]
        hs = [l.s01/l.s11 for l in A.loci]
        M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m, ones(L))
        P,_= fixedpointit(M, ones(L));
        pm = P[end,:,1];
        ps = singlelocuseq(M);
        ds = abs.(ps .- pm)
        push!(sss, ss...)
        push!(hss, hs...)
        push!(wss, pm...)
        push!(dss, ds...)
    end
    sss, hss, wss, dss
end 


# Exponential-uniform model
# -------------------------
Ns = 20
s  = 0.01
Ne = Ns/s
k  = 5
N  = _Ne2N(Ne, k)
u  = s * 0.005
n  = 10
Lss = [0.5, 1., 1.5]
mms = [0.1, 0.3, 0.5]
ss  = 0.005:0.005:0.08
hs  = 0:0.05:1

out = map(Lss) do Ls
    L = ceil(Int, Ls / s)
    map(mms) do ms
        E1, Y1 = msedfe(dfe, ss, hs, L, N, k, u, ms*s, n=n)
        @info Ls, ms
        E1, Y1
    end
end

ED = mapreduce(x->first.(x), hcat, out)
pp = mapreduce(x->last.(x), hcat, out)

pp = map(zip(Lss,out)) do (Ls,x)
    map(zip(mms, last.(x))) do (m,y)
        ttl = "\$L\\bar{s}=$Ls, m/\\bar{s}=$m\$"
        heatmap(ss .- minimum(ss), hs, y, colorbar=false, title=ttl, titlefont=9,
                xlabel="\$s\$", ylabel="\$h\$")
    end 
end |> x->vcat(x...)
plot(pp..., size=(650,600))


Ns = 20
s  = 0.01
Ne = Ns/s
k  = 5
N  = _Ne2N(Ne, k)
u  = s * 0.005
κ = 1
λ = 1/s
Lss = [0.5, 1., 1.5]
mms = 0.05:0.15:0.5
dfe1 = IndependentDFE(Gamma(κ, 1/λ), Beta(1, 1))

out11 = tmap(Lss) do Ls
    L = ceil(Int, Ls / s)
    n = 150000/L
    map(mms) do ms
        @info Ls, ms, L, n
        ss, hs, ws, df = msedfe1(dfe1, L, N, k, u, ms*s, n=n)
    end
end

MP = marginalplot(out11, Lss, mms, dfe1, yh=(0,1.3))

function jointplot(out, Lss, mms, xmx=0.1)
    ps = map(zip(out, Lss)) do (X, Ls)
        map(zip(X, mms)) do ((ss, hs, ws, ds), m)
            D = @sprintf "%.2f" mean(ws)
            ttl = "\$L\\bar{s}=$Ls, m/\\bar{s} = $m, D = $D\$"
            histogram2d(-ss, hs, weights=Weights(ws), colorbar=false,
                        bins=(0:0.005:xmx,20), title=ttl, titlefont=7)
        end 
    end |> x->hcat(x...) 
    ps = permutedims(ps)
    plot(ps..., size=(600,580), xlim=(0,xmx*1.01), tickfont=7, xlabel="\$s\$", ylabel="\$h\$")
end

function marginalplot(out, Lss, mms, dfe, xmx=0.1)
    ps = map(zip(out, Lss)) do (X, Ls)
        map(zip(X, mms)) do ((ss, hs, ws, ds), m)
            D = @sprintf "%.2f" mean(ws)
            ttl = "\$L\\bar{s}=$Ls, m/\\bar{s} = $m, D = $D\$"
            p1 = stephist(-ss, weights=Weights(ws), bins=0:0.005:xmx,
                           title=ttl, titlefont=7, color=:gray, fill=true, alpha=0.5,
                           normalize=true, xlabel="\$s\$")
            plot!(0:0.005:xmx, s->pdf(dfe.sd, s), color=:black)
            p2 = stephist(hs, weights=Weights(ws), bins=0:0.05:1,
                           color=:gray, fill=true, alpha=0.5,
                           normalize=true, xlabel="\$h\$")
            plot!(0:0.01:1, h->MultilocusIsland.pdfh(dfe, h), color=:black)
            plot(p1, p2, legend=false)
        end 
    end |> x->hcat(x...) 
    ps = permutedims(ps)
    plot(ps..., size=(900,380), tickfont=5)
end

function marginalplot(out, Lss, mms, dfe; xmx=0.1, yh=(0,Inf), kwargs...)
    ps = map(zip(out, Lss)) do (X, Ls)
        P1 = plot(0:0.001:xmx, s->pdf(dfe.sd, s), color=:black, xlim=(0,xmx), label="")
        P2 = plot(0:0.01:1, h->MultilocusIsland.pdfh(dfe, h), color=:black, xlim=(0,1), ylim=yh)
        map(enumerate(zip(X, mms))) do (i,((ss, hs, ws, ds), m))
            Hs = [-hs ; hs ; 2 .- hs]
            Ws = [ws  ; ws ; ws]
            kd = kde(Hs, weights=Weights(Ws))
            D = @sprintf "%.2f" mean(ws)
            ttl = "\$L\\bar{s}=$Ls\$"
            lab = @sprintf "%.2f" m
            p1 = density!(P1, -ss, weights=Weights(ws), normalize=true,
                          title=ttl,
                          label="\$m/\\bar{s} = $lab\$", xlabel="\$s\$", color=i,
                          legend=:topright; kwargs...)
            p2 = plot!(P2, 0:0.01:1, h->pdf(kd, h)*3, legend=false,
                          xlabel="\$h\$", color=i; kwargs...)
        end 
        plot(P1, P2, layout=(2,1))
    end |> x->hcat(x...) 
    ps = permutedims(ps)
    plot(ps..., size=(550,280), tickfont=7, layout=(1,3))
end
marginalplot(out11, Lss, mms, dfe1, trim=true)

function scatterplot(out, Lss, mms, xmx=0.1)
    ps = map(zip(out, Lss)) do (X, Ls)
        map(zip(X, mms)) do ((ss, hs, ws, ds), m)
            D = @sprintf "%.2f" mean(ws)
            ttl = "\$L\\bar{s}=$Ls, m/\\bar{s} = $m, D = $D\$"
            za = get.(Ref(cs), ws)
            p1 = scatter(-ss, hs, color=za, alpha=0.2, ms=1, title=ttl, titlefont=7,
                         xlabel="\$s\$", ylabel="\$h\$")
            zb = get.(Ref(cs), hs)
            p2 = scatter(-ss, ws, color=zb, alpha=0.2, ms=1, xlabel="\$s\$", ylabel="\$\\mathbb{E}[p]\$")
            plot(p1, p2, legend=false)
        end 
    end |> x->hcat(x...) 
    ps = permutedims(ps)
    plot(ps..., size=(900,380), tickfont=5)
end
scatterplot(out11, Lss, mms)

#serialize("data/dfes.jls", [(dfe1, out11), (dfe2, out21), (dfe3, out31)])

out = deserialize("data/dfes.jls")

map(enumerate(out)) do (i, (d, o))
    JP = jointplot(o, Lss, mms)
    savefig("$pth/joint$i.svg")
    MP = marginalplot(o, Lss, mms, d)
    savefig("$pth/marginal$i.svg")
    SP = scatterplot(o, Lss, mms)
    savefig("$pth/scatter$i.png")
end

Ps = map(enumerate(out)) do (i, (d, o))
end

MP = marginalplot(out[1][2], Lss, mms, out[1][1], yh=(0,1.3))
MP = marginalplot(out[2][2], Lss, mms, out[2][1], yh=(0,4.5))
MP = marginalplot(out[3][2], Lss, mms, out[3][1], yh=(0,3.5))


# good plot?
dfe1, out11 = out[1]
Ps = map(zip(Lss, out11)) do (Ls, x)
    map(zip(mms, x)) do (m, y)
        ss, hs, ws, ds = y    
        idx = rand(1:length(ws), 1000)
        scatter(ws[idx] .- ds[idx], ds[idx], color=get.(Ref(cs), hs[idx]))
    end
end
Ps = vcat(Ps...)
plot(Ps...)



# Logistic model
# --------------

dfe2 = Logisticsbyh(Gamma(κ, 1/λ), (s/2, 0.5), (0.1, 0.99), 1.)
out21 = tmap(Lss) do Ls
    L = ceil(Int, Ls / s)
    map(mms) do ms
        @info Ls, ms
        ss, hs, ws, df = msedfe1(dfe2, L, N, k, u, ms*s, n=1000)
    end
end

PP2 = doplot1(out21, Lss, mms)

# Caballero & Keightley model
# ---------------------------

out3 = map(Lss) do Ls
    L = ceil(Int, Ls / s)
    map(mms) do ms
        E1, Y1 = msedfe(dfe, ss, hs, L, N, k, u, ms*s, n=n)
        @info Ls, ms
        E1, Y1
    end
end

dfe3 = CKGamma(κ, λ, 1/3)
out31= map(Lss) do Ls
    L = ceil(Int, Ls / s)
    map(mms) do ms
        msedfe1(dfe3, L, N, k, u, ms*s, n=1000)
    end
end

doplot1(out31, Lss, mms)


# Differentiation compared to single-locus
# ---------------------------------------

s̄ = 0.01
κ = 1.
λ = κ/s̄
dfe1 = IndependentDFE(Gamma(κ, 1/λ), Beta(1, 1))
u   = 0.005*s̄
nrep= 2
mms = 0.05:0.05:0.7
Lss = 0.4:0.1:2
Nss = [4, 8, 16]
k   = 5
Xs = map(Nss) do Ns
    N  = _Ne2N(Ns/s̄, k) 
    map(mms) do m
        Ps = map(Lss) do Ls
            L  = ceil(Int, Ls/s̄)
            @info L, m
            ds = tmap(1:nrep) do _
                A  = Architecture([randlocus(dfe1) for i=1:L])
                sa = [l.s01 for l in A.loci]
                sb = [l.s11 - 2l.s01 for l in A.loci]
                hs = [l.s01/l.s11 for l in A.loci]
                M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
                P,_= fixedpointit(M, ones(L));
                pm = P[end,:,1];
                aps = singlelocuseq(M);
                D = sum(pm .- ps)
            end
            mean(ds)
        end
        vcat(Ps...)
    end |> x->hcat(x...)
end 

Z = hcat([collect(Lss ./ s̄) for i=1:size(Xs[1], 2)]...)
Ys = map(X->X ./ Z, Xs)
maximum.(Ys)

Ps = map(zip(Ys, Nss)) do (X, Ns)
    heatmap(mms, Lss, X, size=(300,260), colorbar=false, 
            title="\$N_e\\bar{s} = $Ns\$", clim=(0,0.4),
            xlabel="\$m/\\bar{s}\$", ylabel="\$L\\bar{s}\$")
end
cb = heatmap((0:0.01:1).*ones(101,1), colorbar=false, xticks=false,
             yticks=(1:25:101, string.(0:0.1:0.4)), framestyle=:default)
Py = plot(Ps..., cb, layout=grid(1,4, widths=[0.33,0.33,0.33,0.01]), size=(700,190), 
     bottom_margin=3Plots.mm, left_margin=2Plots.mm, tickfont=7)
# savefig("$pth/driftdiff.svg")


# Variance and gff
# ================
L = 100
s̄ = 0.01
nrep = 50
Nes = 10
k = 5
N = _Ne2N(Nes/s̄, k)
u = s̄*0.005
ms = [0.05, 0.1, 0.2]
κs = [8, 4, 2, 1, 1/2, 1/4]
res = map(ms) do m
    ps = map(κs) do κ
        @info κ
        map(1:nrep) do _
            df = IndependentDFE(Gamma(κ, s̄/κ), Dirac(0.5)) #Beta(1,1))
            A  = Architecture([randlocus(df) for i=1:L])
            hs = [l.s01/l.s11 for l in A.loci]
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
            P,_= fixedpointit(M, ones(L));
            pm = P[end,:,1];
            ps = singlelocuseq(M);
            mean(pm), mean(ps)
        end 
    end
    A  = Architecture(DipLocus(-s̄/2, -s̄), L)
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
    P,_= fixedpointit(M, ones(1));
    pm = P[end,1,1];
    A  = Architecture(DipLocus(-s̄/2, -s̄), 1)
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
    P,_= fixedpointit(M, ones(1));
    p1 = P[end,1,1];
    pm, p1, ps
end

P3 = plot()
m = 0.1
ps = map(κs) do κ
    df = IndependentDFE(Gamma(κ, s̄/κ), Dirac(0.5)) #Beta(1,1))
    A  = Architecture([randlocus(df) for i=1:L])
    smn, smx = extrema([-l.s11 for l in A.loci])
    smn = @sprintf "%.3f" smn
    smx = @sprintf "%.3f" smx
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
    P,_= fixedpointit(M, ones(L));
    pm = P[end,:,1];
    plot!(sort(pm, rev=true), line=:steppost, legend=:topright, 
          label="\$\\kappa=$κ\$", ylabel="\$\\mathbb{E}[p]\$", xlabel="locus")
end
A  = Architecture(DipLocus(-s̄/2, -s̄), L)
M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
P,_= fixedpointit(M, ones(1));
pm = P[end,1,1];
A  = Architecture(DipLocus(-s̄/2, -s̄), 1)
M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
P,_= fixedpointit(M, ones(1));
p1 = P[end,1,1];
hline!([pm], ls=:solid, color=:black, alpha=0.5)
hline!([p1], ls=:dash, color=:black, alpha=0.5)
P3 = plot(P3, title="(C)", legend=false, tickfont=7)

P2 = plot(title="(B)")
map(κ->plot!(Gamma(κ, s̄/κ), label="\$\\kappa = $κ\$"), κs)
plot!(ylabel="density", tickfont=7, xlim=(0,0.05), ylim=(0,200), legend=:topright, xlabel="\$s\$", legendfont=6)

al = ["8", "4", "2", "1", "\$\\frac{1}{2}\$", "\$\\frac{1}{4}\$"]
P1 = plot(title="(A)")
n1 = length(κs)
n2 = length(ms)
zs = map(enumerate(zip(ms, res))) do (j,(m, xys))
    y1, y2, xs = xys
    annotate!((j-1)*n1 + 0.6, 1.04, text("\$m/\\bar{s}=$m\$", 9, :left))
    ys = map(enumerate(zip(κs, xs))) do (i,(α, x_))
        x = first.(x_)
        y = last.(x_)
        # boxplot of mean barrier-wide 'multilocus predicted' differentiation across replicates
        boxplot!([6*(j-1)+i], x, label="", fillalpha=0.9, lw=1, ms=2, whisker_width=0, color=i)
        # boxplot of mean barrier-wide 'single locus predicted' differentiation across replicates
        #violin!([6*(j-1)+i], y, label="", fillalpha=0.5, lw=0, ms=2, outliers=false, whisker_width=0, color=i)
        plot!([0.5+6*(j-1)+i-1,0.5+6*(j-1)+i], [mean(y),mean(y)], color=i)
    end
end
map(enumerate(zip(ms, res, zs))) do (j,(m, xys, z))
    y1, y2, xs = xys
    plot!((0.5+6*(j-1)):1:(6*j+0.5), fill(y1,7), color=:black, ls=:solid, alpha=0.5)
    plot!((0.5+6*(j-1)):1:(6*j+0.5), fill(y2,7), color=:black, ls=:dash, alpha=0.5)
end
vline!((n1+0.52):n1:n2*n1, color=:black, ylim=(0,1), label="")
plot!(xticks=(1:n1*n2, repeat(al, n2)), xtickfont=7, xlim=(0.5,n1*n2 + 0.5),
      xlabel="\$\\kappa\$", ylabel="\$\\mathbb{E}[p]\$", size=(350,250),
      legend=false, top_margin=3Plots.mm)

plot(P1, plot(P2, P3, layout=(2,1)), titlefont=7,
     layout=grid(1,2,widths=[0.7,0.3]), size=(520,260))



# Var[h] variation
# ----------------
L = 100
s̄ = 0.01
nrep = 50
Nes = 10
k = 5
N = _Ne2N(Nes/s̄, k)
u = s̄*0.005
ms = [0.05, 0.1, 0.2, 0.3]
αs = [8, 4, 2, 1, 1/2, 1/4]
res = map(ms) do m
    ps = map(αs) do α
        @info α
        map(1:nrep) do _
            df = IndependentDFE(Dirac(s̄), Beta(α,α))
            A  = Architecture([randlocus(df) for i=1:L])
            hs = [l.s01/l.s11 for l in A.loci]
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
            P,_= fixedpointit(M, ones(L));
            pm = P[end,:,1];
            mean(pm)
        end 
    end
    A  = Architecture(DipLocus(-s̄/2, -s̄), L)
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
    P,_= fixedpointit(M, ones(1));
    pm = P[end,1,1];
    A  = Architecture(DipLocus(-s̄/2, -s̄), 1)
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
    P,_= fixedpointit(M, ones(1));
    p1 = P[end,1,1];
    pm, p1, ps
end

al = ["8", "4", "2", "1", "\$\\frac{1}{2}\$", "\$\\frac{1}{4}\$"]
plot()
n1 = length(αs)
n2 = length(ms)
map(enumerate(zip(ms, res))) do (j,(m, xys))
    y1, y2, xs = xys
    annotate!((j-1)*n1 + 0.6, 1.04, text("\$m/\\bar{s}=$m\$", 9, :left))
    map(enumerate(zip(αs, xs))) do (i,(α, x))
        boxplot!(x, label="", fillalpha=0.9, lw=1, ms=2, whisker_width=0)
    end
end
map(enumerate(zip(ms, res))) do (j,(m, xys))
    y1, y2, xs = xys
    plot!((0.5+6*(j-1)):1:(6*j+0.5), fill(y1,7), color=:black, ls=:solid, alpha=0.5)
    plot!((0.5+6*(j-1)):1:(6*j+0.5), fill(y2,7), color=:black, ls=:dash, alpha=0.5)
end
vline!((n1+0.52):n1:n2*n1, color=:black, ylim=(0,1), label="")
plot!(xticks=(1:n1*n2, repeat(al, n2)), xtickfont=7, xlim=(0.5,n1*n2 + 0.5),
      xlabel="\$\\alpha\$", ylabel="\$\\mathbb{E}[p]\$", size=(350,250),
      legend=false, top_margin=3Plots.mm)

# Swamping thresholds
L = 100
s̄ = 0.01
Nes = 10
k = 5
N = _Ne2N(Nes/s̄, k)
u = s̄*0.005
ms = 0:0.02:0.9
κs = [8, 4, 2, 1, 1/2, 1/4]
ps = tmap(κs) do κ
    df = IndependentDFE(Gamma(κ, s̄/κ), Beta(1,1))
    A  = Architecture([randlocus(df) for i=1:L])
    res = map(ms) do m
        M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
        P,_= fixedpointit(M, ones(L));
        pm = P[end,:,1];
    end
    @info κ
    X = hcat(res...)
end

plot(legend=:topright, size=(300,220))
map(zip(κs,ps)) do (κ,X)
    plot!(ms, vec(sum(X, dims=1)), ylim=(0, L), label="\$\\kappa = $κ\$")
end
plot!(xlim=(0,0.9), ylabel="\$\\mathbb{E}[D]\$", xlabel="\$m/\\bar{s}\$")

function eqgff(p, pq, A)
    map(1:length(p)) do i
        @unpack sa, sb = MultilocusIsland.sasb(A[i])
        exp(2*(sa*p[i] + sb*pq[i]))
    end
end

P4 = plot()
m = 0.05
ps = map(κs) do κ
    df = IndependentDFE(Gamma(κ, s̄/κ), Beta(1,1))
    A  = Architecture([randlocus(df) for i=1:L])
    smn, smx = extrema([-l.s11 for l in A.loci])
    smn = @sprintf "%.3f" smn
    smx = @sprintf "%.3f" smx
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
    P,_= fixedpointit(M, ones(L));
    p  = P[end,:,1];
    pq = P[end,:,2]
    gg = eqgff(p, pq, A)
    g = @sprintf "%.2f" prod(gg)
    plot!(sort(gg), label="\$g=$g, \\kappa=$κ\$", ylabel="\$g\$",
          xlabel="locus", legend=:bottomright, size=(300,220))
end
plot(P4)



