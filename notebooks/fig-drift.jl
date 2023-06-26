using MultilocusIsland, StatsBase, ColorSchemes, ThreadTools, Optim
using Plots, PlotThemes, Distributions, Parameters; theme(:hokusai)
cs = ColorSchemes.viridis
cs = ColorSchemes.Hokusai3

# 1. Nes by Ls for given m/s
# ==========================
s   = 0.02
k   = 50
Lss = 0.04:0.04:2
Nes = 1:1:20
hs  = [0., 0.5, 1.]
mms = [0.1, 0.2, 0.4, 0.8]

res =
    map(mms) do ms
        map(hs) do h
            map(Nes) do Nes_
                map(Lss) do Ls
                    L  = ceil(Int, Ls/s)
                    Ne = Nes_ / s
                    N  = _Ne2N(Ne, k)
                    A  = Architecture(DipLocus(-s*h, -s), L)
                    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, Ne=Ne, u=s*0.005, A=A), ms*s, ones(L))
                    P,_= fixedpointit(M, [1.])
                    p  = P[end,1,1]
                    pq = P[end,1,2]
                    ℓ = MultilocusIsland.eqload(M, p, pq)
                    ℓmx = 1-exp(-Ls)
                    ℓ/ℓmx
                    p
                end
            end |> x->hcat(x...)
        end
    end

Ps = map(enumerate(hs)) do (j,h)
    P = heatmap(1:4length(Nes), Lss, hcat(getindex.(res, j)...), 
            xticks=(5:5:4length(Nes), repeat(5:5:maximum(Nes), 4)))
    vline!(20.5:20:80, color=:white, legend=false, clim=(0,1))
    for (i,k)=enumerate(5:20:4length(Nes))
        annotate!(P, k-4, 2.15, text("\$m/s=$(mms[i])\$", 8, :left))
    end
    annotate!(P, 70, 0.3, text("\$h=$h\$", 8, :left, :white))
    annotate!(P, 94, 1.0, text("\$\\mathbb{E}[p]\$", 8, :center, :black))
    plot(P, size=(600,100), tickfont=6, ylabel="\$Ls\$", 
         xlabel=j==3 ? "\$N_es\$" : "", 
         right_margin=4Plots.mm)
end
plot(Ps..., layout=(3,1), size=(450,380))


# 2. m/s by Ls for various Nes
# ============================
Lss = 0.04:0.04:2
Nes = [2, 4, 8, 16]
hs  = [0., 0.5, 1.]
mms = 0.05:0.05:1
res =
    map(Nes) do Nes_
        map(hs) do h
            map(mms) do ms
                @show ms, h, Nes_
                map(Lss) do Ls
                    L  = ceil(Int, Ls/s)
                    Ne = Nes_ / s
                    N  = _Ne2N(Ne, k)
                    A  = Architecture(DipLocus(-s*h, -s), L)
                    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, Ne=Ne, u=s*0.005, A=A), ms*s, ones(L))
                    P,_= fixedpointit(M, [1.])
                    p  = P[end,1,1]
                    pq = P[end,1,2]
                    ℓ = MultilocusIsland.eqload(M, p, pq)
                    ℓmx = 1-exp(-Ls)
                    ℓ/ℓmx
                    p
                end
            end |> x->hcat(x...)
        end
    end

Ps = map(enumerate(hs)) do (j,h)
    K = length(mms)
    P = heatmap(1:4K, Lss, hcat(getindex.(res, j)...), 
            xticks=(5:5:4K, repeat(0.25:0.25:maximum(mms), 4)))
    vline!(20.5:20:80, color=:white, legend=false, clim=(0,1))
    for (i,k)=enumerate(5:20:4K)
        annotate!(P, k-4, 2.18, text("\$N_es=$(Nes[i])\$", 8, :left))
    end
    annotate!(P, 70, 0.3, text("\$h=$h\$", 8, :left, :white))
    annotate!(P, 94, 1.0, text("\$\\mathbb{E}[p]\$", 8, :center, :black))
    plot(P, size=(600,100), tickfont=6, ylabel="\$Ls\$", 
         xlabel=j==3 ? "\$m/s\$" : "", 
         right_margin=4Plots.mm)
end
plot(Ps..., layout=(3,1), size=(450,380))


# 3. m/s by Nes for various Ls
# ============================
Lss = [0.5, 1., 1.5, 2.]
Nes = 1:1:20
hs  = [0., 0.5, 1.]
mms = 0.05:0.05:1
res =
    map(Lss) do Ls
        map(hs) do h
            map(mms) do ms
                @show ms, h, Ls
                map(Nes) do Nes_
                    L  = ceil(Int, Ls/s)
                    Ne = Nes_ / s
                    N  = _Ne2N(Ne, k)
                    A  = Architecture(DipLocus(-s*h, -s), L)
                    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, Ne=Ne, u=s*0.005, A=A), ms*s, ones(L))
                    P,_= fixedpointit(M, [1.])
                    p  = P[end,1,1]
                    pq = P[end,1,2]
                    ℓ = MultilocusIsland.eqload(M, p, pq)
                    ℓmx = 1-exp(-Ls)
                    ℓ/ℓmx
                    p
                end
            end |> x->hcat(x...)
        end
    end

Ps = map(enumerate(hs)) do (j,h)
    K = length(mms)
    P = heatmap(1:4K, Nes, hcat(getindex.(res, j)...), 
                xticks=([0.5 ; 5:5:4K], [0. ; repeat(0.25:0.25:maximum(mms), 4)]))
    vline!(20.5:20:80, color=:white, legend=false, clim=(0,1))
    for (i,k)=enumerate(5:20:4K)
        annotate!(P, k-4, 23, text("\$Ls=$(Lss[i])\$", 8, :left))
    end
    annotate!(P, 71, 2.5, text("\$h=$h\$", 8, :left, :white))
    annotate!(P, 96, 10.5, text("\$\\mathbb{E}[p]\$", 8, :center, :black))
    plot(P, size=(600,100), tickfont=6, ylabel="\$N_es\$", 
         xlabel=j==3 ? "\$m/s\$" : "", 
         right_margin=4Plots.mm)
end
PH = plot(Ps..., layout=(3,1), size=(450,380), top_margin=2Plots.mm)


# 4. Allele frequency distributions to show
# =========================================
s   = 0.02
k   = 5 
Ls  = 1.0
Nes = [2, 4, 8, 16]
hs  = [0., 0.5, 1.]
ms  = 0.2

res2 =
    map(enumerate(hs)) do (i,h)
        tmap(enumerate(Nes)) do (j,Nes_)
            L  = ceil(Int, Ls/s)
            Ne = Nes_ / s
            N  = _Ne2N(Ne, k)
            A  = Architecture(DipLocus(-s*h, -s), L)
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, Ne=Ne, u=s*0.005, A=A), ms*s, ones(L))
            P,_= fixedpointit(M, [1.])
            p  = P[end,1,1]
            pq = P[end,1,2]
            xs = expectedsfs(M, p, pq, step=0.025, f=log10)
            _,Q= simulate(M, 60000, zeros(L), drop=10000, thin=10)
            @show ms, h, Ls
            xs, Q
        end
    end

pls = map(zip(hs, res2)) do (h, X)
    pl = plot(title="\$h=$h\$", titlefont=8,
              xlabel=h==1.0 ? "\$p\$" : "", 
              ylabel="\$\\log_{10}\\phi(p)\$",
              legend=h==0.0 ? :bottomright : false)
    map(enumerate(zip(Nes, X))) do (k,(Ne, Y))
        xs, Q = Y
        x, y = xs[1]
        plot!(x, y, color=k, label="\$N_es = $(Nes[k])\$", legendfont=6)
        x, y = sfs(1 .- vec(Q), step=0.05, f=log10)
        scatter!(x, y, color=k, ms=2, label="")
    end
    pl
end
annotate!(pls[1], -0.49, 2.8, text("\$\\mathrm{(A)}\$", 9, :left))
annotate!(pls[1],  1.12, 2.8, text("\$\\mathrm{(B)}\$", 9, :left))
PD = plot(pls..., layout=(3,1), tickfont=6, ylim=(-5,1.3))

plot(PD, PH, layout=grid(1,2,widths=[0.20,0.8]), size=(640,380))

# 5. Drift vs deterministic illustrations
# =======================================
L   = 150
Ls  = 1.5
Nes = [1.,2.,4.,8.,16.]
k   = 5
s   = Ls/L
hs  = [-1., 0., 0.25, 0.75, 1, 2]
xs  = vec(collect(Iterators.product(Nes, hs)))
ms2 = 0.05:0.02:1.25


Z = map(xs) do (Nes, h)
    Ne  = Nes / s
    N = _Ne2N(Ne, k)
    A = Architecture(DipLocus(-s*h, -s), L)
    y2 = map(ms2) do ms
        M = MainlandIslandModel(HapDipDeme(Ne=Ne, N=N, k=k, u=s*0.005, A=A), ms*s, ones(L))
        P,_= fixedpointit(M, [1.])
        pm = P[end,1,1]
        dpm = solve(M)
        ms, dpm, pm
    end
    Nes, h, y2
end

dd = Dict(h=>plot() for h=hs)
map(1:length(Z)) do i
    (N,h,y2) = Z[i]
    yl = (h ∈ hs[[1,4]]) ? "\$\\mathbb{E}[p]\$" : ""
    xl = h > hs[3] ? "\$m/s\$" : ""
    xs, ys = first.(y2), last.(y2)
    plot!(dd[h], xs, ys, lw=1, xlabel=xl, ylabel=yl, title="\$h=$h\$",
          legend=h==hs[end] ? :bottomright : false, label="\$N_es=$N\$")
    N == 16 ? plot!(xs, getindex.(y2, 2), color=:black, label="\$N_es \\rightarrow \\infty\$") : nothing
end
plot([dd[h] for h in hs]..., size=(600,300))


# 6. Drift in haplodiplontics
# ===========================

L   = 50
Ls  = 1.0
se  = Ls/L
Nes = 5
ks  = [1, 5, 10] 
s   = Ls/L
u   = s*0.005
he  = 0.9
τs  = [0.2, 0.8]
mms = sort([0.05:0.05:0.35; [0.175, 0.225]]) 

res3 = tmap(ks) do k
    map(τs) do τ
        s = se/(2-τ)
        h = 1 - (1-he*(2-τ))/τ
        N = _Ne2N(Nes/se, k)
        A = Architecture(HapDipLocus(-s*(1-τ), -s*h*τ, -s*τ), L)
        map(mms) do ms
            @info k, N, τ, ms
            M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), ms*se, ones(L))
            _,Q= simulate(M, 150000, zeros(L), drop=50000, thin=10)
            (k, N, τ, ms, Q)
        end
    end
end

res4 = begin
    A = Architecture(DipLocus(-se*he, -se), L)
    mmss = 0:0.005:maximum(mms)
    ys = map(mmss) do ms
        M = MainlandIslandModel(HapDipDeme(Ne=Nes/se, N=0, k=0, u=u, A=A), ms*se, ones(L))
        P,_= fixedpointit(M, [1.])
        P[end,1,1]
    end
    mmss, ys
end

P = plot(title="\$L=$L, N_es_e=$Nes, h_e=$he, u/s=0.005\$")
map(zip(ks,res3)) do (k, X)
    ys = map(zip(τs, X)) do (t, Y)
        qs = 1 .- map(mean, last.(Y))
        N  = Y[1][2] 
        scatter!(mms, qs, label="\$k=$k, N=$N, \\tau=$t\$")
    end
end
plot!(res4, label="diffusion", color=:black, alpha=0.5)
plot(P, xlabel="\$m/s_e\$", ylabel="\$\\mathbb{E}[p]\$", legend=:bottomleft,
     size=(330,250), legendfont=7)

# --------------------------------------------------------------------------------------------
Nes = [2,4,8,16]
PP = map([0, 0.5, 1.0]) do h
    Pl = plot(title="\$h=$h\$", ylabel= h==0 ? "\$\\log_{10}\\phi(p)\$" : "")
    tmap(enumerate(Nes)) do (i,Ns)
        Ne  = Ns / s
        N = _Ne2N(Ne, k)
        A = Architecture(DipLocus(-s*h, -s), L)
        ms = 0.2
        M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=s*0.005, A=A), ms*s, ones(L))
        _,Q=simulate(M, 60000, zeros(L), drop=10000, thin=10)
        P,_= fixedpointit(M, [1.])
        p = P[end,:,1]
        pq = P[end,:,2]
        xs_, ys_ = sfs(vec(Q), step=0.05, f=log10)
        ys = expectedsfs(M, p, pq, step=0.005, f=log10)
        plot!(Pl, ys, label="\$N_es = $Ns\$", color=i, legend= h==0 ? :bottomleft : false)
        scatter!(Pl, xs_, reverse(ys_), label="", color=i, markerstrokecolor=i, ms=3)
    end
    Pl
end

plot(PP..., ms=2, layout=(1,3), xlabel="\$p\$", size=(630,170),
     left_margin=3Plots.mm, bottom_margin=4Plots.mm, legendfont=7)
        

# Detailed look at effects of drift with dominance
Ns = 8.
h  = 0.25
Ne  = Ns / s
N = _Ne2N(Ne, k)
A = Architecture(DipLocus(-s*h, -s), L)
ms = 0.2
M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=s*0.005, A=A), ms*s, ones(L))
_,Q=simulate(M, 51000, zeros(L), drop=1000, thin=5)
P,_= fixedpointit(M, [1.])
p = P[end,:,1]
pq = P[end,:,2]
xs_, ys_ = sfs(vec(Q), step=0.025, f=log10)
ys = expectedsfs(M, p, pq, step=0.005, f=log10)
plot(ys, label="\$N_es = $Ns\$", color=1, legend= h==0 ? :bottomleft : false)
scatter!(xs_, reverse(ys_), label="", color=1, markerstrokecolor=1, ms=3)


# Obtain the migration rate for which some critcial 𝔼p is obtained.
using Optim, Printf
L   = 50
Ls  = 1.2
s   = Ls/L
Nes = 1:1:32
k   = 500
hs  = -0.25:0.25:1.25
pc  = 0.05

# objective function
function obj(pc, Nes, s, h, L, k=5)
    Ne = Nes / s
    N  = _Ne2N(Ne, k)
    A  = Architecture(DipLocus(-s*h, -s), L)
    function ff(m)
        M  = MainlandIslandModel(HapDipDeme(N=N, k=k, Ne=Ne, u=s*0.005, A=A), m*s, ones(L))
        P,_= fixedpointit(M, [1.])
        abs(P[end,1,1] - pc)
    end
end

#hs = [0.6, 0.5, 0.4]
hs = [0., 0.5, 1.]
Nes = 1:0.2:12
xs = map(hs) do h
    tmap(Nes) do Ne
        res = optimize(obj(pc, Ne, s, h, L), 0., 3.)
        (Ne, h, res.minimizer)
    end |> x->vcat(x...) 
end

P1 = plot(title="\$Ls=$Ls, L=$L\$")
map(enumerate(xs)) do (i,x)
    h = x[1][2]
    Nes = first.(x)
    xss = Nes 
    c = get(cs, (h - minimum(hs)) / (maximum(hs) - minimum(hs)))
    plot!(P1, xss, last.(x), label="\$h=$(@sprintf "%.2f" h)\$", marker=false,
          ms=2, color=c, markerstrokecolor=c)
#    hline!([1-h], color=c, ls=:dot, label="")
end
plot(P1, xlabel="\$N_e s\$", ylabel="\$m_c/s \\ [p_c=$pc]\$",
     size=(300,270), legendfont=8, legend=:topleft, ylim=(0,1.))

# closer look at what happens at the critical Ne...
ms = 0.01:0.02:1.
hs  = [0., 0.5, 1]
plot()
Ness = [1, 2, 4, 8, 16, 32, 64]
map(hs) do h
    PP = plot(legend = h == 0 ? :topright : false, 
              legendfont=7,
              xlabel="\$m/s\$", ylabel="\$\\mathbb{E}[p]\$", title="\$h=$h\$")
    map(Ness) do Nes
        Ne = Nes / s
        N  = _Ne2N(Ne, k)
        ps = map(ms) do m
            A  = Architecture(DipLocus(-s*h, -s), L)
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, Ne=Ne, u=s*0.005, A=A), m*s, ones(L))
            P,_= fixedpointit(M, [1.])
            P[end,1,1]
        end
        cc = log2.(Ness)
        c = get(cs, (log2(Nes) - minimum(cc))/maximum(cc))
        plot!(ms, ps, color=c, label="\$N_es = $Nes\$")
    end
    hline!([pc], color=:gray, label="")
end |> x->plot(x..., layout=(1,3), size=(720,200), left_margin=3Plots.mm, 
               bottom_margin=3Plots.mm)

ms = 0.01:0.01:1.
hs  = [1]
plot()
Ness = 1:32
map(hs) do h
    PP = plot(legend = h == 0 ? :topright : false, 
              legendfont=7,
              xlabel="\$m/s\$", ylabel="\$\\mathbb{E}[p]\$", title="\$h=$h\$")
    map(Ness) do Nes
        Ne = Nes / s
        N  = _Ne2N(Ne, k)
        ps = map(ms) do m
            A  = Architecture(DipLocus(-s*h, -s), L)
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, Ne=Ne, u=s*0.005, A=A), m*s, ones(L))
            P,_= fixedpointit(M, [1.])
            P[end,1,1]
        end
        cc = log2.(Ness)
        c = get(cs, (log2(Nes) - minimum(cc))/maximum(cc))
        plot!(ms, ps, color=c, label="\$N_es = $Nes\$")
    end
    hline!([pc], color=:gray, label="")
end  |> first |> x->plot(x, size=(350,250))



# This doesn't appear to be interesting: determining pc from deterministic
# theory, and matching it with stochastic.
hs  = [0., 0.25, 0.5, 0.75, 1]
Nes = [1, 2, 4, 8, 16, 32, 64, 128]
γs = [0.25, 0.5, 0.75, 0.9]
pls = map(γs) do γ
    ys = map(hs) do h
        map(Nes) do Ne
            s = Ls/L
            sa = -s*h
            sb = -s + 2s*h 
            x, ps, _ = MultilocusIsland.critical_ms(sa, sb, L, tol=1e-9)
            @show h, Ne, ps[1]
            pc = γ*ps[1]
            res = optimize(obj(pc, Ne, s, h, L), 0., 3.)
            (Ne, h, res.minimizer)
        end |> x->vcat(x...) 
    end
    plot()
    map(enumerate(ys)) do (i,x)
        h = x[1][2]
        c = get(cs, h)
        plot!(first.(x), last.(x), label="\$h=$(@sprintf "%.2f" h)\$", marker=true,
              ms=3, color=c, markerstrokecolor=c)
    end
    plot!(xlabel="\$N_e s\$", 
          ylabel="\$m_c/s\$",
          title="\$\\gamma = $γ\$",
          size=(300,230), legend=:bottomright)
end 
plot(pls..., layout=(1,4), size=(900,220), margin=3Plots.mm)


# ------------
# Migration rate for which we arrive at γ * deterministic critical p.
hs  = [0., 0.25, 0.5, 0.75, 1]
Nes = [1, 2, 4, 8, 16, 32, 64, 128]
γ   = 0.9
Ls  = 1.2
L   = 50
s   = Ls/L
k   = 5
u   = s*0.005

function getcriticalm(s, h, L, γ)
    sa = -s*h
    sb = -s + 2s*h 
    x, ps, _ = MultilocusIsland.critical_ms(sa, sb, L, tol=1e-9)
    pc = γ*ps[1]
    ps[1], pc
end

function obj(pc, Nes, s, h, L, k)
    Ne = Nes / s
    N  = _Ne2N(Ne, k)
    A  = Architecture(DipLocus(-s*h, -s), L)
    function ff(m)
        M  = MainlandIslandModel(HapDipDeme(N=N, k=k, Ne=Ne, u=s*0.005, A=A), m*s, ones(L))
        P,_= fixedpointit(M, [1.])
        abs(P[end,1,1] - pc)
    end
end

ys = map(hs) do h
    pcd, pc = getcriticalm(Ls/L, h, L, γ)
    mcs = map(Ns->optimize(obj(pc, Ns, Ls/L, h, L, k), 0., 3.).minimizer, Nes)
    pcd, pc, mcs
end

P1 = plot()
map(enumerate(ys)) do (i,x)
    h = hs[i]
    c = get(cs, h)
    plot!(Nes, x[end], label="\$h=$(@sprintf "%.2f" h)\$", marker=true,
          ms=3, color=c, markerstrokecolor=c)
end
plot!(xlabel="\$N_e s\$", 
      ylabel="\$m_c/s\$",
      title="\$\\gamma = $γ\$",
      size=(300,230), legend=:bottomright)

ms = 0.01:0.01:1.2
pps = map(enumerate(hs)) do (i,h)
    plot()
    pcd, pc, mcs = ys[i]
    map(Nes) do Nes_
        Ne = Nes_ / s
        N  = _Ne2N(Ne, k)
        ps = map(ms) do m
            A  = Architecture(DipLocus(-s*h, -s), L)
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, Ne=Ne, u=s*0.005, A=A), m*s, ones(L))
            P,_= fixedpointit(M, [1.])
            P[end,1,1]
        end
        plot!(ms, ps, color=i, label="")
    end
    scatter!([mcs], [pc], color=i, label="", markerstrokecolor=i, legend=false,
              xlabel="\$m/s\$", ylabel="\$\\mathbb{E}[p]\$")
end  

P2 = plot(P1, pps..., layout=(2,3), size=(800,400))

# second derivative?    
function swampfirstdiff(ps, ms)
    yy = map(i->(ps[i+1] - ps[i])/step(ms), 1:length(ps)-1)
    zz = map(i->(yy[i+1] - yy[i])/step(ms), 1:length(yy)-1)
    start = 1
    while (zz[start] > 0) && (start < length(yy)-1)
        start += 1
    end
    idx = all(zz .> 0) ? 1 : start + argmin(yy[start:end]) - 1
    return ms[idx], ps[idx]
end

Xs = map(enumerate(hs)) do (i, h)
    PP = plot()
    xs = map(Nes) do Nes_
        Ne = Nes_ / s
        N  = _Ne2N(Ne, k)
        ps = map(ms) do m
            A  = Architecture(DipLocus(-s*h, -s), L)
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, Ne=Ne, u=u, A=A), m*s, ones(L))
            P,_= fixedpointit(M, [1.])
            P[end,1,1]
        end
        mc, pc = swampfirstdiff(ps, ms)
        plot!(ms, ps, color=i, label="")
        scatter!([mc], [pc], color=i, label="", markerstrokecolor=i, legend=false,
              xlabel="\$m/s\$", ylabel="\$\\mathbb{E}[p]\$")
        (Nes_, mc, pc)
    end
    plot!(ms, u ./ ((ms .* s) .+ u), color=:gray)
    xs, PP
end 
plot(last.(Xs)...)

P2 = plot()
for (i, xs) in enumerate(first.(Xs))
    plot!(first.(xs), getindex.(xs, 2), marker=true, markerstrokecolor=i, label="\$h=$(@sprintf "%.2f" hs[i])\$")
end
P2 = plot(P2, xlabel="\$N_es\$", ylabel="\$m_c/s\$", legend=:bottomright)

plot(P2, last.(Xs)..., legendfont=7, size=(700,340))

# -----------------
L   = 40
Ls  = 0.8
Nes = [2.,4.,8.,16.,32]
k   = 5
s   = Ls/L
hs  = [0., 0.5, 1.]
xs  = vec(collect(Iterators.product(Nes, hs)))
ms1 = 0.01:0.01:0.85

map(hs) do h
    @show h
    A = Architecture(DipLocus(-s*h, -s), L)
    PP = plot()
    map(enumerate(Nes)) do (i,Ns)
        xs = map(ms1) do ms
            Ne = Ns/s
            N = _Ne2N(Ne, k)
            MM = MainlandIslandModel(HapDipDeme(N=N, k=k, u=s*0.005, A=A), ms*s, ones(L))
            P,_= fixedpointit(MM, [1.])
            Epm = P[end,1,1]
            pm = solve(MM)
            (ms, pm, Epm)
        end
        plot!(getindex.(xs, 1), getindex.(xs, 3), color=i, 
              markerstrokecolor=i, 
              label="\$N_es=$Ns\$")
        i == length(Nes) && plot!(getindex.(xs, 1), getindex.(xs, 2), color=:black, 
              markerstrokecolor=i, 
              label="\$N_es=\\infty\$")
    end
    PP
end |> x->plot(x...)

Z = map(xs) do (Nes, h)
    Ne  = Nes / s
    N = _Ne2N(Ne, k)
    A = Architecture(DipLocus(-s*h, -s), L)
    y2 = map(ms1) do ms
        M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=s*0.005, A=A), ms*s, ones(L))
        P,_= fixedpointit(M, [1.])
        pm = P[end,1,1]
        ms, pm
    end
    Nes, h, y2
end

dd = Dict(h=>plot() for h=hs)
map(1:length(Z)) do i
    (N,h,y2) = Z[i]
    yl = (h ∈ hs[[1,4]]) ? "\$\\mathbb{E}[p]\$" : ""
    xl = h > hs[3] ? "\$m/s\$" : ""
    xs, ys = first.(y2), last.(y2)
    plot!(dd[h], xs, ys, lw=2, xlabel=xl, ylabel=yl, title="\$h=$h\$",
          legend=h==hs[end] ? :bottomright : false, label="\$N_es=$N\$")
end
plot([dd[h] for h in hs]..., size=(600,300))




#---------------------
# Fix Ls?
Ls = 1.2
L  = [10:1:25 ; 25:5:500]
ms = [0.001, 0.005, 0.01]
Ns = [100, 200, 400, 800, 1600, 3200]
h  = 0.5

lss = [:solid, :dash, :dot]
plot()
map(enumerate(ms)) do (j,m)
    xx = map(Ns) do N
        map(L) do L
            s = Ls/L
            @show N, s
            A  = Architecture(DipLocus(-s*h, -s), L)
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=s*0.005, A=A), m, ones(L))
            P,_= fixedpointit(M, [1.])
            P[end,1,1]
        end
    end |> x->hcat(x...)
    map(enumerate(Ns)) do (i,N)
        plot!(Ls ./ L, xx[:,i], label="\$N=$N\$", xlabel="\$s\$",
              ylabel="\$\\mathbb{E}[p]\$", color=i, ls=lss[j])
    end
end
plot!()

Ls = 2.5
L  = 50
s  = Ls/L 
u  = s*0.005
A0 = Architecture(DipLocus(-s/2, -s), L)
ss = 10 .^ -range(3, 1, 100)
h  = 1.
m  = s/4
M  = MainlandIslandModel(HapDipDeme(N=1000, k=k, u=u, A=A0), m, ones(length(A0)))
P,_= fixedpointit(M, [1.]); P

PP = map([0, 0.5, 1.]) do h
    X = map(Ns) do N
        map(ss) do s
            A = vcat(A0, Architecture(DipLocus(-s*h, -s), 1))
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m, ones(length(A)))
            P,_= fixedpointit(M, [1., 1.])
            P[end,end,1]
        end 
    end|> x->hcat(x...)
    plot(ylim=(0,1))
    map(enumerate(Ns)) do (i,N) 
        plot!(ss, X[:,i], label="\$N=$N\$")
    end
    #plot!(ss, max.(0, 1 .- 2m ./ ss), color=:black, label="single locus", alpha=0.4)
    plot!(xscale=:log10, title="\$h=$h, Ls=$Ls, L=$L\$")
end 

plot(PP..., layout=(1,3), size=(700,200), legend=:topleft,
     xlabel="\$s_{L+1}\$", ylabel="\$\\mathbb{E}[p]\$",
     margin=3Plots.mm)


# The gff
# =======
L   = 50
k   = 5
Lss = [0.5, 1., 1.5]
hs  = 0:0.2:1
Nss = [2, 4, 8, 16]

Ps = map(Lss) do Ls
    s = Ls/L
    u = s*0.005
    ms = 0.01:0.01:0.8
    PP = map(Nss) do Ns
        PP = plot(title="\$Ls = $Ls, N_es = $Ns\$", legend= Ls==Lss[1] && Ns == Nss[1] ? :topright : false)
        map(hs) do h
            Ne = Ns/s
            N  = _Ne2N(Ne, k)
            A  = Architecture(DipLocus(-s*h, -s), L)
            sa = -s*h
            sb = -s*(1 - 2h)
            gs = map(ms) do m
                M  = MainlandIslandModel(HapDipDeme(N=N, k=k, Ne=Ne, u=u, A=A), m*s, ones(L))
                P,_= fixedpointit(M, [1.])
                Ep = P[end,1,1]
                Epq = P[end,1,2]
                g = exp(2L*(sa*Ep + sb*Epq))
            end
            plot!(ms, gs, label="\$h=$h\$")
        end
        PP
    end
end |> x->vcat(x...)
plot(Ps..., layout=(3,4), xlabel="\$m/s\$", ylabel="\$g\$")


# Quantify drift load vs. migration load?
# =======================================
L   = 50
k   = 5
Ls  = 1.1
h   = 0.4
s   = Ls/L
u   = s*0.005
mss = 0.0:0.01:0.8
Ns  = 10.
N   = _Ne2N(Ns/s, k)
A   = Architecture(DipLocus(-s*h, -s), L)

ys  = map(mss) do ms
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, Ne=Ns/s, u=u, A=A), ms*s, ones(L))
    P,_= fixedpointit(M, [1.])
    p1, pq1 = P[end,:,:]
    l1 = MultilocusIsland.eqload(M, p1, pq1)
    p2 = solve(M)
    l2 = 1-MultilocusIsland.meanfitness(A[1], p2)^L
    p1, p2, l1, l2, (l1 - l2)/l1
end

plot(mss, first.(ys))
plot!(mss, getindex.(ys, 2))
plot!(mss, getindex.(ys, 3))
plot!(mss, getindex.(ys, 4))
plot!(twinx(), mss, getindex.(ys,5), color=:black, ylim=(0,1))
hline!([1-exp(-Ls)])

Lss = 0.5:0.05:1.5
Nss = 1:32
xs = map(0.05:0.1:0.55) do ms
    map(Lss) do Ls
        s = Ls/L
        u = s*0.005
        map(Nss) do Ns
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, Ne=Ns/s, u=u, A=A), ms*s, ones(L))
            P,_= fixedpointit(M, [1.])
            p1, pq1 = P[end,:,:]
            l1 = MultilocusIsland.eqload(M, p1, pq1)
            p2 = solve(M)
            l2 = 1-MultilocusIsland.meanfitness(A[1], p2)^L
            #p1, p2, l1, l2
            (l1 - l2)/l1
        end 
    end |> x->hcat(x...)
end

plot(heatmap.(Ref(Lss), Ref(Nss), xs)..., size=(900,370))


