using MultilocusIsland
using Plots, PlotThemes, Printf, ColorSchemes; theme(:hokusai)
cs = ColorSchemes.viridis

# could we work this out a bit more carefully? Like, when does LD matter as a
# function f(Ls, sf)...?

# a focal locus in an additive background
Lss= 0.6:0.2:1.6
L  = 50
u  = 0.0001
Ns = 20.
Ne = 1000 #Ns/s
k  = 5
N  = _Ne2N(Ne, k) 
sf = 0.05

P1s = map([1, 0]) do h
    mss = 0:0.005:0.65
    Af = Architecture(DipLocus(-h*sf, -sf), 1)
    PP = plot()
    map(enumerate(Lss)) do (i,Ls)
        s = Ls/L
        Ab = Architecture(DipLocus(-0.5s, -s ), L)
        A = vcat(Ab, Af)
        ps = map(mss) do ms
            M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), ms*sf)
            P,_ = fixedpointit(M, ones(2))
            P[end,:,1]
        end |> x->hcat(x...) |> permutedims 
        sl = @sprintf "%.2f" s
        lab = "\$Ls=$Ls\$"
        P1 = plot!(mss, ps, color=[i i], alpha=[0.3 1], lw=[2 1.5],
                  label=["" lab], legend=:topright,
                  xlabel="\$m/s_f\$", ylabel="\$\\mathbb{E}[p]\$", 
                  title="\$h=$h\$")
    end
    ps2 = map(mss) do ms
        M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=Af), ms*sf)
        P,_= fixedpointit(M, ones(1));
        P[end,1,1]
    end
    plot!(PP, mss, ps2, label="single locus", color=:black)
    vline!([1/6, 1/4], color=:black, ls=:dash, alpha=0.4, label="")
    PP 
end
plot(P1s..., size=(600,200))

P2s = map([1, 0]) do h
    Af = Architecture(DipLocus(-h*sf, -sf), 1)
    map(zip([6,4], [:solid, :solid])) do (sm, ls)
        PP = plot(ylim=(-15,2))
        m = sf/sm
        M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=Af), m)
        P,_= fixedpointit(M, ones(1));
        p  = P[end,:,1]
        pq = P[end,:,2]
        ys2 = expectedsfs(M, p, pq, step=0.005, f=log10)
        plot!(ys2, color=:gray, alpha=0.4, lw=2, label="single locus", xlabel="\$p\$",
              ylabel=(sm == 4) ? "\$\\log_{10}\\phi(p)\$" : "", ls=ls,
              size=(300,250), legend=false)
        map(enumerate(Lss)) do (i,Ls)
            s = Ls/L
            Ab = Architecture(DipLocus(-0.5s, -s ), L)
            A = vcat(Ab, Af)
            M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m)
            P,_ = fixedpointit(M, ones(2))
            P[end,:,1]
            p  = P[end,:,1]
            pq = P[end,:,2]
            ys = expectedsfs(M, p, pq, step=0.005, f=log10)
            plot!(PP, ys[2], color=i, lw=2, label="\$Ls=$Ls\$", title="\$m=s_f/$sm\$", ls=ls)
        end
        PP
    end
end |> x->vcat(x...)
plot(P2s..., size=(600,200), layout=(1,4), ylim=(-20, 2))

P11 = plot(P1s[1], P2s[1:2]..., size=(500,200), layout=grid(1,3, widths=[0.4,0.3,0.3]), margin=1Plots.mm)
P12 = plot(P1s[2], P2s[3:4]..., size=(500,200), layout=grid(1,3, widths=[0.4,0.3,0.3]), margin=1Plots.mm, legend=false)
PP  = plot(P11, P12, layout=(2,1), size=(630,340), legendfont=7, titlefont=9)
for (x, p) in zip("ABCDEF", PP.subplots)
    p.attr[:title] = "($x) " * p.attr[:title]
end
plot(PP)
   
