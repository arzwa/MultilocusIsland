using MultilocusIsland
using Plots, PlotThemes, Printf, ColorSchemes; theme(:hokusai)
cs = ColorSchemes.viridis

function critical_point(h, Ls)
    if h == 0.5
        ms = exp(Ls - 1)/(2Ls)
        pc = 1 - 1/Ls
    elseif h == 0.
        ms = (0.25 + √(Ls*(Ls + 8))/(4Ls))^2 * exp(Ls/4 + √(Ls*(Ls+8))/4 -1)
        pc = (3Ls - √(Ls*(Ls + 8)))/4Ls
    elseif h == 1.
        pc = 0.5 + Ls/4
        ms = (0.25 - Ls^2/16)*exp((Ls*(Ls^2 + 4Ls + 4))/8)
    else
        @warn "not implemented" 
        pc = NaN
        ms = NaN
    end
    if !(0 < pc < 1)
        pc, ms = NaN, NaN
    end
    ms, pc
end

# Dominance in the multilocus model: deterministic case
PP = map([0., 0.5, 1.0]) do h
    P = plot()
    map(enumerate([0.01,0.25,0.5,0.75,1.,1.25,1.5])) do (i,Ls)
        s  = 0.01
        sa = -s*h
        sb = -s -2sa
        x, y, x_, y_, = findroots_ms(sa, sb, Ls/s)
        color = i == 1 ? :lightgray : i-1
        plot!(x,y,   color=color, lw=3, alpha=0.7, label="\$Ls=$(@sprintf "%.2f" Ls)\$")
        plot!(x_,y_, color=color, alpha=0.7, label="", ls=:dot)
        ms, pc = critical_point(h, Ls)
        cc = :firebrick
        scatter!([ms], [pc], color=cc, label="", markerstrokecolor=cc)
    end
    P1 = plot(P, title="\$s=0.01, h=$(@sprintf "%.2f" h)\$", xlabel="\$m/s\$", xlim=(0,1.25),
              ylabel="\$\\tilde{p}\$", legend=h > 0.5 ? :bottomright : :topright)
end

hs = [-1,0,0.5,1,2]
hs = 0:0.25:1
xs = let s=0.01
    #map(0:0.2:1) do h
    map(hs) do h
        sa = -s*h
        sb = -s + 2s*h
        xs = []
        ys = []
        ms0 = 0.
        for Ls=0.01:0.01:1.5
            x, ps, _ = MultilocusIsland.critical_ms(sa, sb, Ls/s, tol=1e-9)
            push!(xs, (x, ps[1]))
            push!(ys, critical_point(h, Ls)[1])
        end
        xs, ys
    end
end

P2 = plot(xlabel="\$Ls\$", ylabel="\$(m/s)_c\$")
map(enumerate(zip(hs, xs))) do (i,(h,x))
    plot!(0.01:0.01:1.5, first.(x[1]), label="\$h=$(@sprintf "%.2f" h)\$", lw=2, color=i)
    plot!(0.01:0.01:1.5, x[2], label="", lw=2, color=i, ls=:dot)
end
P2 = plot!(legend=:topleft, ylim=(0,1.25), legendfont=7, bg_legend=:white)

plot(PP..., P2, size=(550,450))
#savefig("$pth/detdom.svg")
        

# genetic load (assume stable equilibrium)
function load(p, s11, sb, L)
    q = 1 - p
    1 - exp(L*(s11*q - sb*p*q))
end

msmax = 1.25
css = [1,3,5]
PP = map(enumerate([0.01,0.75,1.25])) do (i,Ls)
    P = plot(title=Ls==0.01 ? "single locus" : "\$Ls=$(@sprintf "%.2f" Ls)\$", 
             legend=Ls==1.5 ? :topleft : :bottomright)
    map(zip([0., 0.5, 1.0], css)) do (h, col)
        s  = 0.01
        L  = Ls/s 
        sa = -s*h
        sb = -s -2sa
        x, p, _, _, = findroots_ms(sa, sb, Ls/s)
        x = [x; [x[end] + 0.0001, msmax]]
        p = [p; [0., 0.]]
        y = load.(p, -s, sb, L) 
        plot!(x, y, color=col, lw=3, alpha=0.7, 
              label="\$h = $(@sprintf "%.2f" h)\$")
    end
    hline!([1-exp(-Ls)], ls=:dot, color=:black, label="\$1-e^{-Ls}\$")
    P
end
plot(PP..., ylabel="load", xlabel="\$m/s\$", size=(700,180), layout=(1,3), margin=3Plots.mm)


# Equilibrium
eqf(p, Ls, h, ms, q=1-p) = q*(h*p + (1-2h)*p*q) - exp(-2Ls*(h*p + (1-2h)*p*q))*ms*p 

map([0, 0.5, 1]) do h
    map([0.5, 1.0, 1.5]) do Ls
        P = plot()
        for (i,ms) in enumerate([0.01, 0.05, 0.1, 0.2, 0.5])
            plot!(map(p->(p, eqf(p, Ls, h, ms)), 0:0.01:1), color=i,
                  title="\$Ls=$Ls, h=$h\$", label="\$m/s=$ms\$", 
                  legend=Ls==1.5 && h==0.0 ? :bottomleft : false)
        end
        P
    end |> x->plot(x..., layout=(1,3))
end |> x->plot(x..., layout=(3,1), size=(750,600), framestyle=:origin)


eqf2(p, Ls, h, ms, q=1-p) = h*q + (1-2*h)*q^2 - ms*exp(-2Ls*(h*p + (1-2h)*p*q))

map([0, 0.5, 1]) do h
    map([0.5, 1.0, 1.5, 3.0]) do Ls
        P = plot()
        for (i,ms) in enumerate([0.1, 0.2, 0.5, 1.5, 3.0])
            plot!(map(p->(p, eqf2(p, Ls, h, ms)), 0:0.01:1), color=i,
                  title="\$Ls=$Ls, h=$h\$", label="\$m/s=$ms\$", 
                  legend=Ls==3.0 && h==0.0 ? :topright : false)
        end
        P
    end |> x->plot(x..., layout=(1,4), ylim=(-1,Inf))
end |> x->plot(x..., layout=(3,1), size=(820,600), framestyle=:origin)


f3(p, Ls, h, ms, q=1-p) = h*q + (1-2*h)*q^2 
me(p, Ls, h, ms, q=1-p) = ms*exp(-2Ls*(h*p + (1-2h)*p*q))
map([0, 0.5, 1]) do h
    map([0.5, 1.0, 1.5, 3.0]) do Ls
        P = plot()
        for (i,ms) in enumerate([0.1, 0.2, 0.5, 1.5, 3.0])
            plot!(map(p->(p, f3(p, Ls, h, ms)), 0:0.01:1), color=i,
                  title="\$Ls=$Ls, h=$h\$", label="\$m/s=$ms\$", 
                  legend=Ls==3.0 && h==0.0 ? :topright : false)
            plot!(map(p->(p, me(p, Ls, h, ms)), 0:0.01:1), color=i,
                  ls = :dot)
        end
        P
    end |> x->plot(x..., layout=(1,4), ylim=(-1,Inf))
end |> x->plot(x..., layout=(3,1), size=(820,600), framestyle=:origin)
