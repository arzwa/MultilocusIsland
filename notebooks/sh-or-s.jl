using MultilocusIsland, Plots, Distributions, PlotThemes; theme(:hokusai)

# For a single locus

function equilibria(m, s, h)
    sa = -s*h
    sb = -s*(1-2h)
    sb == 0 && return 1+m/sa, NaN
    4sb/m > (sa/m)^2 && return NaN, NaN
    D  = √((sa/m)^2 - 4sb/m)
    q1 = (-sa + m*D)/(2sb)
    q2 = (-sa - m*D)/(2sb)
    dq1 = sa - m + 2*(sb - sa)*q1 - 3sb*q1^2
    dq2 = sa - m + 2*(sb - sa)*q2 - 3sb*q2^2
    return (1-q2, 1-q1)
end

P1 = plot()
m = 0.1
mss = 0:0.001:1
map(enumerate(0.1:0.2:0.9)) do (i,h)
    ys = map(mss) do msh  #m/sh
        sh = m / msh
        s  = sh / h
        equilibria(m, s, h)
    end
    plot!(mss, first.(ys), lw=3, color=i, label="\$h=$h\$")
    plot!(mss, last.(ys), ylim=(0,1), color=i, alpha=0.2, lw=2, label="")
end
#plot!(mss, x->1-x, color=:black, ls=:dot)
plot!(xlabel="\$m/sh\$", ylabel="\$\\tilde{p}\$")

P2 = plot()
m = 0.1
map(enumerate(0.1:0.2:0.9)) do (i,h)
    ys = map(mss) do ms  #m/s
        s  = m/ms
        equilibria(m, s, h)
    end
    plot!(mss, first.(ys), lw=3, color=i, label="\$h=$h\$")
    plot!(mss, last.(ys), ylim=(0,1), color=i, alpha=0.2, lw=2, label="")
end
plot!(xlabel="\$m/s\$", ylabel="\$\\tilde{p}\$")

plot(P1, P2, size=(600,200))


P1 = plot()
m = 0.1
mss = 0:0.0001:0.1
map(enumerate(0.1:0.2:0.9)) do (i,h)
    ys = map(mss) do msh  #m/sh
        sh = m / msh
        s  = sh / h
        equilibria(m, s, h)
    end
    plot!(mss, 1 .- first.(ys), lw=3, color=i, label="\$h=$h\$")
    plot!(mss, 1 .- last.(ys), ylim=(0,0.12), color=i, alpha=0.2, lw=2, label="")
end
#plot!(mss, x->1-x, color=:black, ls=:dot)
plot!(xlabel="\$m/sh\$", ylabel="\$\\tilde{q}\$")

P2 = plot()
m = 0.1
map(enumerate(0.1:0.2:0.9)) do (i,h)
    ys = map(mss) do ms  #m/s
        s  = m/ms
        equilibria(m, s, h)
    end
    plot!(mss, 1 .- first.(ys), lw=3, color=i, label="\$h=$h\$")
    plot!(mss, 1 .- last.(ys), ylim=(0,0.3), color=i, alpha=0.2, lw=2, label="")
end
plot!(xlabel="\$m/s\$", ylabel="\$\\tilde{q}\$")

plot(P1, P2, size=(600,200), margin=3Plots.mm)
