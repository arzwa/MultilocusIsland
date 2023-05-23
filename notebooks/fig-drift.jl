using MultilocusIsland, StatsBase, ColorSchemes, ThreadTools
using Plots, PlotThemes, Distributions, Parameters; theme(:hokusai)
cs = ColorSchemes.viridis
cs = ColorSchemes.Hokusai3

L   = 40
Ls  = 0.8
Nes = [1.,2.,4.,8.,16.]
k   = 5
s   = Ls/L
hs  = [-1., 0., 0.25, 0.75, 1, 2]
xs  = vec(collect(Iterators.product(Nes, hs)))
ms2 = 0.05:0.1:0.55


Z = map(xs) do (Nes, h)
    Ne  = Nes / s
    N = _Ne2N(Ne, k)
    A = Architecture(DipLocus(-s*h, -s), L)
    y2 = map(ms1) do ms
        M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=s*0.005, A=A), ms*s, ones(L))
        P,_= fixedpointit(M, [1.])
        pm = P[end,1,1]
        dpm = solve(MM)
        ms, pm, dpm
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
xs = map(hs) do h
    tmap(Nes) do Ne
        res = optimize(obj(pc, Ne, s, h, L), 0., 3.)
        (Ne, h, res.minimizer)
    end |> x->vcat(x...) 
end

P1 = plot(title="\$Ls=$Ls, L=$L\$")
map(enumerate(xs)) do (i,x)
    h = x[1][2]
    c = get(cs, (h - minimum(hs)) / (maximum(hs) - minimum(hs)))
    plot!(P1, first.(x), last.(x), label="\$h=$(@sprintf "%.2f" h)\$", marker=true,
          ms=2, color=c, markerstrokecolor=c)
#    hline!([1-h], color=c, ls=:dot, label="")
end
plot(P1, xlabel="\$N_e s\$", ylabel="\$m_c/s \\ [p_c=$pc]\$",
     size=(450,330), legendfont=6, ylim=(0,1.5))


# This doesn't appear to be interesting: determining pc from deterministic
# theory, and matching it with stochastic.
ys = map(hs) do h
    map(Nes) do Ne
        s = Ls/L
        sa = -s*h
        sb = -s + 2s*h 
        x, ps, _ = MultilocusIsland.critical_ms(sa, sb, L, tol=1e-9)
        @show h, Ne, ps[1]
        res = optimize(obj(ps[1], Ne, s, h, L), 0., 3.)
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
plot!(xlabel="\$N_e s\$", ylabel="\$m_c/s \\ [p_c=$pc]\$", size=(300,230), legend=:bottomright)


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
