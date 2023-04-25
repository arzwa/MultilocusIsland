using MultilocusIsland
using Plots, PlotThemes, Printf, ColorSchemes; theme(:hokusai)
cs = ColorSchemes.viridis

# single-locus
P = plot()
map(enumerate([-1., 0., 0.5, 1.0, 2.0])) do (i,h)
    s  = 0.01
    sa = -s*h
    sb = -s -2sa
    x, y, x_, y_, = findroots_ms(sa, sb, 1, ms0=0.000001)
    plot!(x,y,   color=i, lw=3, alpha=0.7, label="\$h=$h\$", legend=:topright)
    plot!(x_,y_, color=i, alpha=0.7, label="", ls=:dot)
end
plot(P, size=(300,230), xlabel="\$m/s\$", ylabel="\$\\tilde{p}\$", title="\$w_{Aa} = 1-sh, w_{aa} = 1-s\$")

# diffusion
function ϕ(p, N, m, u, sa, sb) 
    q = 1-p
    p^(2N*u - 1) * q^(2N*u + 2N*m - 1) * exp(N*(2sa*q + sb*q^2)) 
end

s  = 0.01
Ns = 10.
h  = 0.15
u  = 0.005*s
N  = Ns/s
xx = map(0:0.1:0.4) do ms
    pp = plot()
    prange = 0.002:0.01:0.998
    map([0.01, 0.99]) do h
        sa = -s*h
        sb = -s*(1-2h)
        phis = map(p->ϕ(p, N, ms*s, u, sa, sb), prange)
        ys = phis ./ sum(phis)
        plot!(prange, ys, legend=false, framestyle=:default, 
             color=h < 0.5 ? :black : :salmon, yticks=false,
             fill=true, fillalpha=0.2, xticks= ms == 0.4 ? (0:0.2:1) : false, 
             ylabel= ms == 0.2 ? "\$\\phi(p)\$" : "",
             xlabel= ms == 0.4 ? "\$p\$" : "")
        vline!([sum(prange .* ys)], color=h < 0.5 ? :black : :salmon)
    end
    _, ymx = ylims(pp)
    annotate!(0.0, ymx*1.15, text("\$m/s = $ms\$", 7, :left))
end 
plot(xx..., layout=(5,1), size=(200,300), xlim=(0,1), yshowaxis=false)
