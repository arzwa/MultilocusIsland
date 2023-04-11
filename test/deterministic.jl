using MultilocusIsland
using Plots, PlotThemes, Printf, ColorSchemes; theme(:hokusai)
cs = ColorSchemes.viridis

s = 0.02
l = 15
h = [zeros(l); fill(0.5,l); ones(l)]
L = length(h)
τ = 1.
N = 2000
k = 5
u = s*0.005
m = 0.01

# Haplodiplontic model
A = Architecture([HapDipLocus(-s*(1-τ), -s*h[i]*τ, -s*τ) for i=1:L], fill(0.5, L)) 
    
ps = map(0:0.01:1) do ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, arch=A, u=u)
    fixedpointit(M, ones(3))[end,:,1]
end |> x->hcat(x...) |> permutedims 
plot(0:0.01:1, ps)

M = HapDipMainlandIsland(N=N, k=k, m=0.2*s, arch=A, u=u)
Ep = fixedpointit(M, ones(3))[end,:,1]
solve(M, [0.99, 0.99, 0.99])

mss = 0.05:0.01:0.9
ps = map(mss) do ms
    @show ms
    M = HapDipMainlandIsland(N=N, k=k, m=ms*s, arch=A, u=u)
    q = fixedpointit(M, ones(3))[end,:,1]
    [q ; solve(M, q)[1]]
end |> x->hcat(x...) |> permutedims 
plot(mss, ps)


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
    end
    P1 = plot(P, title="\$s=0.01, h=$(@sprintf "%.2f" h)\$", xlabel="\$m/s\$",
              ylabel="\$\\tilde{p}\$", legend=h > 0.5 ? :bottomright : :bottomleft)
end

xs = let s=0.01
    map(0:0.2:1) do h
        sa = -s*h
        sb = -s + 2s*h
        xs = []
        ms0 = 0.
        for Ls=0.01:0.01:1.5
            x, ps, _ = MultilocusIsland.critical_ms(sa, sb, Ls/s, tol=1e-9)
            push!(xs, (x, ps[1]))
        end
        xs
    end
end

P2 = plot(xlabel="\$Ls\$", ylabel="\$(m/s)_c\$")
map(zip(0:0.2:1, xs)) do (h,x)
    plot!(0.01:0.01:1.5, first.(x), label="\$h=$(@sprintf "%.2f" h)\$", lw=2)
end
P2 = plot!(legend=:topleft, ylim=(0,1.25), legendfont=7)

plot(PP..., P2, size=(550,450))
savefig("/home/arthur_z/vimwiki/build/img/2023-04-07/detdom.svg")
        

# Note that Himani does not show the unstable equilibrium in her plot 1B!
x, y, x_, y_, = findroots_ms(-0.02, 0, 40)
plot(2x, y, xlim=(0,1.25)); plot!(2x_, y_)
