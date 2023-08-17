# Variance and gff
# ================
using MultilocusIsland, StatsBase, ColorSchemes, ThreadTools, Printf, QuadGK
using Plots, PlotThemes, Distributions, Parameters; theme(:hokusai)
using LogExpFunctions
cs = ColorSchemes.viridis

function randbarrier(df, L, mx=Inf)
    ss = rand(truncated(df.sd, 0.0, mx), L)
    hs = rand(df.hd, L)
    #ss = ss .* mean(df.sd)/(mean(ss))
    #hs = hs .* mean(df.hd)/(mean(hs))
    [DipLocus(-s*h, -s) for (s,h) in zip(ss,hs)]
end

function singlelocuseq(M)
    map(1:length(M.D.A)) do i
        DD = reconstruct(M.D, A=M.D.A[i:i])
        MM = MainlandIslandModel(DD, M.m1)
        P,_= fixedpointit(MM, ones(1))
        P[end,:,1]
    end |> x->first.(x)
end

L = 100
s̄ = 0.01
nrep = 100
Nes = 10
k = 5
N = _Ne2N(Nes/s̄, k)
u = s̄*0.005
ms = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
κs = [8, 4, 2, 1, 1/2, 1/4]
al = ["8", "4", "2", "1", "\$\\frac{1}{2}\$", "\$\\frac{1}{4}\$"]
res = map(ms) do m
    ps = map(κs) do κ
        @info κ
        map(1:nrep) do _
            df = IndependentDFE(Gamma(κ, s̄/κ), Dirac(0.5)) #Beta(1,1))
            #df = IndependentDFE(Gamma(κ, s̄/κ), Beta(1,1))
            A  = Architecture(randbarrier(df, L))
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, 0., ones(L))
            P,_= fixedpointit(M, ones(L));
            pm = P[end,:,1];
            g  = MultilocusIsland.gff(M, pm, P[end,:,2])
            ps = singlelocuseq(M);
            mean(pm), mean(ps), g
        end 
    end
    A  = Architecture(DipLocus(-s̄/2, -s̄), L)
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, 0., ones(L))
    P,_= fixedpointit(M, ones(1));
    pm = P[end,1,1];
    A  = Architecture(DipLocus(-s̄/2, -s̄), 1)
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, 0., ones(L))
    P,_= fixedpointit(M, ones(1));
    p1 = P[end,1,1];
    pm, p1, ps
end

Ps = map(enumerate(zip(ms,res))) do (i,(m, X1))
    P = plot(title="\$m/\\bar{s}=$m\$")
    ms = map(enumerate(zip(κs, X1[3]))) do (j,(k,X2))
        ms = map(x->mean(map(y->y.m ./ (m*s̄), x)), last.(X2))
        boxplot!([j], ms, legend=false, fillalpha=0.9, lw=1, ms=2, whisker_width=0,) 
    end
    P
end
plot(Ps..., xticks=(1:length(al), al), ylim=(0,1),
     xlabel="\$\\kappa\$",ylabel="\$\\bar{g} = \\overline{m_e}/m\$", size=(600,400))

Ps = map(enumerate(collect(zip(ms,res))[2:2:end])) do (i,(m, X1))
    ms = map(enumerate(collect(zip(κs, X1[3]))[1:2:end])) do (j,(k,X2))
        xls = [0.045, 0.085, 0.15]
        P = plot(title="\$m/\\bar{s}=$m\$", xlims=(0,xls[j]))
        colors = [1,3,5]
        map(last.(X2)[1:50]) do xs
            ms = map(y->   y.m ./ (m*s̄), xs)
            ss = map(y-> -2y.sa , xs)
            o = sortperm(ms)
            #scatter!(ss, ms, alpha=0.2, ms=2, color=colors[j])
            plot!(ss[o], ms[o], alpha=0.4, marker=true, ms=2, color=colors[j])
        end
        P
    end
    plot(ms..., ylabel="\$g_i = \\overline{m_{e,i}}/m\$", xlabel="\$s\$", layout=(3,1)) #xlim=(0,0.15)
end
plot(Ps..., legend=false, layout=(1,3), size=(600,600))


# Var[s] = E[s]/κ^2 <=> σ(s) = 0.1 ./ κs = [0.0125, 0.025, 0.05, 0.1, 0.2, 0.4]

P3s = map(zip(["B","C"],[0.1,0.4])) do (ll,m)
    P3 = plot()
    ps = map(κs) do κ
        df = IndependentDFE(Gamma(κ, s̄/κ), Dirac(0.5)) #Beta(1,1))
        A  = Architecture([randlocus(df) for i=1:L])
        ss = [-l.s11 for l in A.loci]
        smn, smx = extrema(ss)
        smn = @sprintf "%.3f" smn
        smx = @sprintf "%.3f" smx
        M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, 0., ones(L))
        P,_= fixedpointit(M, ones(L));
        pm = P[end,:,1];
        oo = sortperm(pm, rev=true)
        plot!(pm[oo], line=:steppost, legend=:topright, 
              label="\$\\kappa=$κ\$", ylabel="\$\\mathbb{E}[p]\$", xlabel="locus")
    end
    A  = Architecture(DipLocus(-s̄/2, -s̄), L)
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, 0., ones(L))
    P,_= fixedpointit(M, ones(1));
    pm = P[end,1,1];
    A  = Architecture(DipLocus(-s̄/2, -s̄), 1)
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, 0., ones(L))
    P,_= fixedpointit(M, ones(1));
    p1 = P[end,1,1];
    hline!([pm], ls=:solid, color=:black, alpha=0.5)
    P3 = plot(P3, title="($ll) \$m/\\bar{s}=$m\$", legend=false, tickfont=7, ylim=(0,1))
end

P2 = plot(title="(B)")
map(κ->plot!(Gamma(κ, s̄/κ), label="\$\\kappa = $κ\$"), κs)
plot!(ylabel="density", tickfont=7, xlim=(0,0.0505), ylim=(0,200), legend=:topright, xlabel="\$s\$", legendfont=6)

al = ["8", "4", "2", "1", "\$\\frac{1}{2}\$", "\$\\frac{1}{4}\$"]
P1 = plot(title="(A)")
n1 = length(κs)
n2 = length(ms)
zs = map(enumerate(zip(ms, res))) do (j,(m, xys))
    y1, y2, xs = xys
    ll = j == 1 ? "\$m/\\bar{s}=$m\$" : "\$$m\$"
    annotate!(j*n1 + 0.6, 1.04, text(ll, 9, :right))
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
end
vline!((n1+0.52):n1:n2*n1, color=:black, ylim=(0,1), label="")
plot!(xticks=(1:n1*n2, repeat(al, n2)), xtickfont=7, xlim=(0.5,n1*n2 + 0.5),
      xlabel="\$\\kappa\$", ylabel="\$\\bar{\\Delta}\$", size=(350,250),
      legend=false, top_margin=3Plots.mm)

PP = plot(P1, plot(P3s..., layout=(2,1)), titlefont=9, margin=2Plots.mm,
     layout=grid(1,2,widths=[0.78,0.22]), size=(650,260))
plot!(PP, inset=(1, bbox(0.36, -0.29, 100Plots.px, 70Plots.px, :center)), subplot=4)
map(κ->plot!(Gamma(κ, s̄/κ), label="\$\\kappa = $κ\$", ylabel="density",
             tickfont=6, xlim=(0,0.0505), guidefont=7, ylim=(0,200), legend=:topright,
             xlabel="\$s\$", legendfont=6, subplot=4), κs)
plot(PP)



# Var[h] variation
# ================
L = 100
s̄ = 0.01
nrep = 50
Nes = 10
k = 5
N = _Ne2N(Nes/s̄, k)
u = s̄*0.005
ms = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
αs = [8, 4, 2, 1, 1/2, 1/4]
res = map(ms) do m
    ps = map(αs) do α
        @info α
        map(1:nrep) do _
            df = IndependentDFE(Dirac(s̄), Beta(α,α))
            A  = Architecture([randlocus(df) for i=1:L])
            hs = [l.s01/l.s11 for l in A.loci]
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, 0., ones(L))
            P,_= fixedpointit(M, ones(L));
            pm = P[end,:,1];
            ps = singlelocuseq(M);
            mean(pm), mean(ps)
        end 
    end
    A  = Architecture(DipLocus(-s̄/2, -s̄), L)
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, 0., ones(L))
    P,_= fixedpointit(M, ones(1));
    pm = P[end,1,1];
    A  = Architecture(DipLocus(-s̄/2, -s̄), 1)
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, 0., ones(L))
    P,_= fixedpointit(M, ones(1));
    p1 = P[end,1,1];
    pm, p1, ps
end

al = ["8", "4", "2", "1", "\$\\frac{1}{2}\$", "\$\\frac{1}{4}\$"]
P1 = plot(title="(A)")
n1 = length(κs)
n2 = length(ms)
zs = map(enumerate(zip(ms, res))) do (j,(m, xys))
    y1, y2, xs = xys
    ll = j == 1 ? "\$m/\\bar{s}=$m\$" : "\$$m\$"
    annotate!(j*n1 + 0.6, 1.04, text(ll, 9, :right))
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
end
vline!((n1+0.52):n1:n2*n1, color=:black, ylim=(0,1), label="")
plot!(xticks=(1:n1*n2, repeat(al, n2)), xtickfont=7, xlim=(0.5,n1*n2 + 0.5),
      xlabel="\$\\alpha\$", ylabel="\$\\bar{\\Delta}\$", size=(350,250),
      legend=false, top_margin=3Plots.mm)


P2 = plot(title="(B)")
map(α->plot!(Beta(α, α), label="\$\\alpha = $α\$"), αs)
plot!(ylabel="density", tickfont=7, xlim=(0,1), ylim=(0,6),
      legend=:outertopright, xlabel="\$h\$", legendfont=6)

P3s = map(zip(["B","C"],[0.1,0.4])) do (ll,m)
    P3 = plot()
    ps = map(αs) do α
        df = IndependentDFE(Dirac(s̄), Beta(α,α))
        A  = Architecture([randlocus(df) for i=1:L])
        ss = [-l.s11 for l in A.loci]
        smn, smx = extrema(ss)
        smn = @sprintf "%.3f" smn
        smx = @sprintf "%.3f" smx
        M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, 0., ones(L))
        P,_= fixedpointit(M, ones(L));
        pm = P[end,:,1];
        oo = sortperm(pm, rev=true)
        plot!(pm[oo], line=:steppost, legend=:topright, 
              label="\$\\alpha=$α\$", ylabel="\$\\mathbb{E}[p]\$", xlabel="locus")
    end
    A  = Architecture(DipLocus(-s̄/2, -s̄), L)
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, 0., ones(L))
    P,_= fixedpointit(M, ones(1));
    pm = P[end,1,1];
    A  = Architecture(DipLocus(-s̄/2, -s̄), 1)
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, 0., ones(L))
    P,_= fixedpointit(M, ones(1));
    p1 = P[end,1,1];
    hline!([pm], ls=:solid, color=:black, alpha=0.5)
    P3 = plot(P3, title="($ll) \$m/\\bar{s}=$m\$", legend=false, tickfont=7, ylim=(0,1))
end

PP = plot(P1, plot(P3s..., layout=(2,1)), titlefont=9, margin=2Plots.mm,
     layout=grid(1,2,widths=[0.78,0.22]), size=(650,260))
plot!(PP, inset=(1, bbox(0.36, -0.29, 100Plots.px, 70Plots.px, :center)), subplot=4)
map(α->plot!(Beta(α, α), label="\$\\alpha = $α\$",
             tickfont=6, xlim=(0,1), guidefont=7, ylim=(0,5), legend=:topright,
             xlabel="\$h\$", legendfont=6, subplot=4), αs)
plot(PP)



