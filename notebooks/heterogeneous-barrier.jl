using MultilocusIsland, StatsBase, ColorSchemes, ThreadTools, Printf
using Plots, PlotThemes, Distributions, Parameters; theme(:hokusai)
cs = ColorSchemes.viridis

# deterministic single-locus predictions
function singlelocuseq_det(M)
    map(1:length(M.D.A)) do i
        @unpack s01, s11 = M.D.A[i]
        sa = s01
        sb = s11 - 2s01
        y, z = MultilocusIsland.roots_stab(sa, sb, M.m, 1)
        if length(y) == 0
            0.
        elseif length(y) == 1
            y[1]
        else
            z[1] < 0 ? y[1] : y[2]
        end
    end
end

# with drift
function singlelocuseq(M)
    map(1:length(M.D.A)) do i
        DD = reconstruct(M.D, A=M.D.A[i:i])
        MM = MainlandIslandModel(DD, M.m, ones(1))
        P,_= fixedpointit(MM, ones(1))
        P[end,:,1]
    end |> x->first.(x)
end

function standardize(xs)
    xmn, xmx = extrema(xs)
    xs .-= xmn
    xs ./ (xmx - xmn)
end

s̄ = 0.01
PPs = map([0.1, 0.3, 0.5]) do m
    Ps = map([0.5,1.0,1.5]) do Ls
        L  = ceil(Int, Ls/s̄)
        κ  = 1.
        λ  = κ/s̄ 
        #dfe= IndependentDFE(Gamma(κ, 1/λ), Beta(1,1)) 
        #dfe = Logisticsbyh(Gamma(κ, 1/λ), (s̄/2, 0.5), (0.1, 0.99), 1.)
        dfe= CKGamma(κ, λ, 1/3)
        k  = 5
        Ns = 10.  
        N  = _Ne2N(Ns/s̄, k) 
        u  = 0.005*s̄
        plot()
        for i=1:10
            A  = Architecture([randlocus(dfe) for i=1:L])
            sa = [l.s01 for l in A.loci]
            sb = [l.s11 - 2l.s01 for l in A.loci]
            hs = [l.s01/l.s11 for l in A.loci]
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
            P,_= fixedpointit(M, ones(L));
            pm = P[end,:,1];
            ps = singlelocuseq(M);
            zs = get.(Ref(cs), hs)
            scatter!(ps, pm, color=zs, markerstrokecolor=zs, ms=3, alpha=0.5)
        end
        plot!(x->x, color=:black, xlim=(0,1), ylim=(0,1),
              titlefont=9,
              xtickfont=7, ytickfont=7,
              title="\$L\\bar{s}=$Ls, m/\\bar{s} = $m\$", 
              xlabel=m  == 0.5 ? "\$\\mathbb{E}[p_i]\$ single locus" : "", 
              ylabel=Ls == 0.5 ? "\$\\mathbb{E}[p_i]\$ multilocus" : "")
    end
    vcat(Ps...)
end 

PP = deepcopy(vcat(PPs...))
plot!(PP[end], repeat([0.9], 20), range(0.1,0.8,20), 
      lw=4, color=get.(Ref(cs), range(0,1,20)))
annotate!(0.82, 0.45, text("\$h\$", 9, :left))
annotate!(0.93, 0.1, text("\$0\$", 9, :left))
annotate!(0.93, 0.8, text("\$1\$", 9, :left))
Py = plot(PP..., legend=false, size=(550,600), alpha=0.2)


#
s̄   = 0.02
nrep= 10
mms = 0.05:0.05:0.7
Lss = 0.4:0.1:2
Nss = [4, 8, 16]
k   = 5
Xs = map(Nss) do Ns
    map(mms) do m
        Ps = map(Lss) do Ls
            L  = ceil(Int, Ls/s̄)
            sd = Exponential(s̄)
            hd = Beta()
            N  = _Ne2N(Ns/s̄, k) 
            u  = 0.005*s̄
            ds = map(1:nrep) do _
                A  = Architecture([randlocus(sd, hd) for i=1:L])
                sa = [l.s01 for l in A.loci]
                sb = [l.s11 - 2l.s01 for l in A.loci]
                hs = [l.s01/l.s11 for l in A.loci]
                M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
                P,_= fixedpointit(M, ones(L));
                pm = P[end,:,1];
                ps = singlelocuseq(M);
                sum(abs.(ps .- pm))
            end
            @info Ls, m
            mean(ds)
        end
        vcat(Ps...)
    end |> x->hcat(x...)
end 

Ps = map(zip(Xs, Nss)) do (X, Ns)
    heatmap(mms, Lss, X ./ L, size=(300,260), colorbar=false, 
            title="\$N_e\\bar{s} = $Ns\$", clim=(0,0.4),
            xlabel="\$m/\\bar{s}\$", ylabel="\$L\\bar{s}\$")
end

cb = heatmap((0:0.01:1).*ones(101,1), colorbar=false, xticks=false,
             yticks=(1:25:101, string.(0:0.1:0.4)), framestyle=:default)
Py = plot(Ps..., cb, layout=grid(1,4, widths=[0.33,0.33,0.33,0.01]), size=(700,190), 
     bottom_margin=3Plots.mm, left_margin=2Plots.mm, tickfont=7)

plot(Py, Px, layout=grid(2,1, heights=[0.2,0.8]), size=(900,900))


# DFE with exponential-uniform input
s̄  = 0.02
Xs = tmap([0.1, 0.3, 0.5, 0.7]) do m
    Ps = map([0.5,1.0,1.5]) do Ls
        L = ceil(Int, Ls/s̄)
        sd = Exponential(s̄)
        hd = Beta()
        k  = 5
        Ns = 50.
        N  = _Ne2N(Ns/s̄, k) 
        u  = 0.005*s̄
        plot()
        sss = Float64[]
        hss = Float64[]
        wss = Float64[]
        for i=1:100
            A  = Architecture([randlocus(sd, hd) for i=1:L])
            sa = [l.s01 for l in A.loci]
            sb = [l.s11 - 2l.s01 for l in A.loci]
            ss = [l.s11 for l in A.loci]
            hs = [l.s01/l.s11 for l in A.loci]
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
            P,_= fixedpointit(M, ones(L));
            pm = P[end,:,1];
            ps = singlelocuseq(M);
            ds = abs.(ps .- pm)
            zs = get.(Ref(cs), hs)
            push!(sss, ss...)
            push!(hss, hs...)
            push!(wss, pm...)
        end
        @show m, Ls
        (m, Ls, sss, hss, wss)
    end 
end 
        
PPs = map(Xs) do X
    map(X) do x
        m, Ls, sss, hss, wss = x
        Ls = ceil(Int, Ls/s̄)
        binsize = s̄/5
        P1 = stephist(-sss, bins=0:binsize:7s̄, weights=Weights(wss), fill=true,
                      normalize=true, color=:lightgray, xlabel="\$s\$")
        plot!(0:0.001:7s̄, map(x->pdf(Exponential(s̄), x), 0:0.001:7s̄), color=:black, 
              titlefont=9,
              xtickfont=6, ytickfont=6,
              title="\$L\\bar{s}=$Ls, m/\\bar{s} = $m\$")
        P2 = stephist(hss, bins=0:0.05:1, weights=Weights(wss), fill=true,
                      normalize=true, color=:lightgray, xlabel="\$h\$")
        plot!(0:0.01:1, map(x->pdf(Beta(1,1), x), 0:0.01:1), color=:black, xtickfont=6, ytickfont=6)
        plot(P1, P2)
    end
end
PP = deepcopy(vcat(PPs...))
Px = plot(PP..., legend=false, size=(900,500), alpha=0.2, layout=(4,3))

PPs = map(Xs) do X
    map(X) do x
        m, Ls, sss, hss, wss = x
        Ls = ceil(Int, Ls/s̄)
        wss ./= sum(wss)
        ρ = pearcor(-sss, hss, wss)
        xmx = round(quantile(-sss, 0.99), digits=2)
        xbn = range(0, xmx, 15)
        ttl = "\$L\\bar{s} = $Ls, m/\\bar{s} = $m\$\n\$\\rho=$(@sprintf "%.2f" ρ)\$"
        P3 = histogram2d(-sss, hss, weights=Weights(wss), bins=(xbn,0:0.06666:1),
                         normalize=:probability, colorbar=false, 
                         xlim=(0,xmx), xticks=range(0,xmx,3), tickfont=7,
                         title=ttl, titlefont=8, 
                         ylims=(0,1),
                         xlabel="\$s\$", ylabel="\$h\$", guidefont=8)
    end
end
Px = plot(vcat(PPs...)..., legend=false, size=(500,780), layout=(4,3), margin=1Plots.mm)

PPs = map(Xs) do X
    map(X) do x
        m, Ls, sss, hss, wss = x
        Ls = ceil(Int, Ls/s̄)
        xmx = round(quantile(-sss, 0.99), digits=2)
        xbn = range(0, xmx, 20)
        ttl = "\$L\\bar{s} = $Ls, m/\\bar{s} = $m\$"
        P1 = scatter(hss, wss, ms=1, color=:gray, alpha=0.2, 
                ylim=(0,1), xlim=(0,1), title=ttl,
                xlabel="\$h\$", ylabel="\$\\mathbb{E}[p]\$")
        P2 = histogram2d(hss, wss, bins=(0:0.05:1, 0:0.05:1), colorbar=false,
                    xlabel="\$h\$", ylabel="\$\\mathbb{E}[p]\$")
        plot(P1, P2)
    end
end
Px = plot(vcat(PPs...)..., legend=false, size=(1000,700), layout=(4,3), margin=1Plots.mm)

PPs = map(Xs) do X
    map(X) do x
        m, Ls, sss, hss, wss = x
        Ls = ceil(Int, Ls/s̄)
        xmx = round(quantile(-sss, 0.99), digits=2)
        xbn = range(0, xmx, 20)
        ttl = "\$L\\bar{s} = $Ls, m/\\bar{s} = $m\$"
        zs = get.(Ref(cs), hss)
        P1 = scatter(-sss, wss, ms=1, color=zs, alpha=0.2, 
                title=ttl,
                xlabel="\$s\$", ylabel="\$\\mathbb{E}[p]\$")
    end
end
Px = plot(vcat(PPs...)..., legend=false, size=(700,700), xscale=:log10,
          layout=(4,3), margin=1Plots.mm, xlim=(5e-5, 1e-1))

PPs = map(Xs) do X
    map(X) do x
        m, Ls, sss, hss, wss = x
        L = ceil(Int, Ls/s̄)
        xmx = round(quantile(-sss, 0.99), digits=2)
        xbn = range(0, xmx, 20)
        ttl = "\$L\\bar{s} = $Ls, m/\\bar{s} = $m\$"
        zs = get.(Ref(cs), wss)
        P1 = scatter(-sss, hss, ms=1, color=zs, alpha=0.5, 
                title=ttl,
                xlabel="\$s\$", ylabel="\$h\$")
    end
end
Py = plot(vcat(PPs...)..., legend=false, size=(500,700), xscale=:log10,
          layout=(4,3), margin=1Plots.mm, xlim=(5e-4, 1e-1))

function pearcor(ss, hs, ws)
    ws ./= sum(ws)
    ms = sum(ss .* ws)
    mh = sum(hs .* ws) 
    cv = sum((ss .- ms) .* (hs .- mh) .* ws)
    vs = sum(((ss .- ms) .^ 2) .* ws)
    vh = sum(((hs .- mh) .^ 2) .* ws) 
    ρ = cv / √(vs*vh) 
end


# dominance rainbow
PPs = map([0.01, 0.05, 0.1, 0.2]) do m
    Ps = map([0.1, 0.5, 1.0, 1.5]) do Ls
        #Ls = 1.0
        s̄  = 0.01
        L  = ceil(Int, Ls/s̄)
        sd = Dirac(s̄)
        hd = Beta()
        k  = 5
        Ns = 8.  
        N  = Ne2N(Ns/s̄, k) 
        u  = 0.005*s̄
        A  = Architecture([DipLocus(-s̄*(i/L), -s̄) for i=1:L])
        hs = [l.s01/l.s11 for l in A.loci]
        M  = HapDipMainlandIsland(N=N, k=k, m=m*s̄, u=u, arch=A)
        xs = summarize_arch(M)
        P,_= fixedpointit(M, ones(xs.K));
        pm = P[end,:,1];
        ps = singlelocuseq(M);
        zs = get.(Ref(cs), hs)
        scatter(ps, pm, color=zs, markerstrokecolor=zs, ms=3)
        plot!(x->x, color=:black, 
              titlefont=8,
              xtickfont=7, ytickfont=7,
              title="\$Ls=$Ls, m/s = $m\$", 
              xlabel=m == 0.2 ? "singlelocus" : "", 
              ylabel=Ls == 0.1 ? "multilocus" : "")
    end
    Ps
end 
plot(vcat(PPs...)..., legend=false, size=(700,700))

        

# predict allele frequency distributions
Ls = 1.0
L  = 50
s̄  = Ls/L
sd = Exponential(s̄)
hd = Beta()
k  = 5
Ns = 10.  
N  = Ne2N(Ns/s̄, k) 
u  = 0.005*s̄
A  = Architecture([randlocus(sd, hd) for i=1:L])
hs = [l.s01/l.s11 for l in A.loci]
m  = 0.20
M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
P,_= fixedpointit(M, ones(L));
p  = P[end,:,1]
pq = P[end,:,2]

ys = expectedsfs(M, p, pq, step=0.005, f=log10)
plot(ys)

_,Q = simulate(M, 210000, zeros(L), drop=10000, thin=10)

PP = plot(size=(450,320))
map(enumerate(1:5)) do (c,i)
    h = A[i].s01/A[i].s11
    s = -A[i].s11
    plot!(ys[i], color=c, label=@sprintf "locus %d \$(s=%.3f, h=%.2f)\$" i s h)
    scatter!(sfs(1 .- Q[:,i], step=0.05, f=log10), markerstrokecolor=c, color=c, label="")
end
PP = plot(PP, legend=:bottomright, ylabel="\$\\log_{10}\\phi(p)\$", xlabel="\$p\$")


mss = 0:0.02:0.8
fps = map(mss) do ms
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), ms*s̄, ones(L))
    P,_= fixedpointit(M, ones(L));
    p  = P[end,:,1]
end
X1 = permutedims(hcat(fps...))

mss2 = 0.05:0.1:0.8
ibs = tmap(mss2) do ms
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), ms*s̄, ones(L))
    _,Q = simulate(M, 210000, zeros(L), drop=10000, thin=10)
    vec(mean(Q, dims=1))
end
X2 = permutedims(hcat(ibs...))

rn = 1:6
PP2 = plot(mss, X1[:,rn], xlabel="\$m/\\bar{s}\$", ylabel="\$\\mathbb{E}[p]\$", legend=false)
scatter!(mss2, 1 .- X2[:,rn], markerstrokecolor=[1 2 3 4 5 6], size=(300,200))
plot!(title="\$L=$L, L\\bar{s} = $Ls, N_e \\bar{s} = $Ns, m/\\bar{s} = $m\$")

plot(PP2, PP, size=(600,230), margin=3Plots.mm, ms=3)




# s by h
L = 100
PPs = map([0.05, 0.1, 0.2]) do m
    Ps = map([0.5,1.0,1.5]) do Ls
        plot()
        for i=1:10
            s̄ = Ls/L
            sd = Exponential(s̄)
            hd = Beta()
            k  = 5
            Ns = 8.  
            N  = Ne2N(Ns/s̄, k) 
            u  = 0.005*s̄
            A  = Architecture([randlocus(sd, hd) for i=1:L])
            sa = [l.s01 for l in A.loci]
            sb = [l.s11 - 2l.s01 for l in A.loci]
            ss = [-l.s11 for l in A.loci]
            hs = [l.s01/l.s11 for l in A.loci]
            M  = HapDipMainlandIsland(N=N, k=k, m=m*s̄, u=u, arch=A)
            xs = summarize_arch(M)
            P,_= fixedpointit(M, ones(xs.K));
            pm = P[end,:,1];
            zs = get.(Ref(cs), pm)
            scatter!(hs, log10.(ss), color=zs, markerstrokecolor=zs, ms=3, alpha=0.5)
        end
        plot!(titlefont=8,
              xtickfont=6, ytickfont=7, ylim=(-5.5,-0.5),
              title="\$L\\bar{s}=$Ls, m/\\bar{s} = $m\$", 
              xlabel=m == 0.2 ? "\$h\$" : "", 
              ylabel=Ls == 0.5 ? "\$\\log_{10}s\$" : "")
    end
    vcat(Ps...)
end 

PP = deepcopy(vcat(PPs...))
plot!(PP[end], range(0.1, 0.9, 20), repeat([-5], 20),  
      lw=4, color=get.(Ref(cs), range(0,1,20)))
annotate!(0.5, -4.7, text("\$\\mathbb{E}[p]\$", 9, :center))
annotate!(0.07, -5, text("\$0\$", 9, :right))
annotate!(0.95, -5, text("\$1\$", 9, :left))
plot(PP..., legend=false, size=(600,650))

L = 100
Xs = map([0.05, 0.1, 0.2]) do m
    Ys = map([0.5,1.0,1.5]) do Ls
        ys = map(1:10) do rep
            s̄ = Ls/L
            sd = Exponential(s̄)
            hd = Beta()
            k  = 5
            Ns = 8.  
            N  = Ne2N(Ns/s̄, k) 
            u  = 0.005*s̄
            A  = Architecture([randlocus(sd, hd) for i=1:L])
            sa = [l.s01 for l in A.loci]
            sb = [l.s11 - 2l.s01 for l in A.loci]
            ss = [-l.s11 for l in A.loci]
            hs = [l.s01/l.s11 for l in A.loci]
            M  = HapDipMainlandIsland(N=N, k=k, m=m*s̄, u=u, arch=A)
            xs = summarize_arch(M)
            P,_= fixedpointit(M, ones(xs.K));
            pm = P[end,:,1];
            collect(zip(pm, ss, hs, sa, sb) )
        end
        (m, Ls, vcat(ys...))
    end
    vcat(Ys...)
end 
Xs = vcat(Xs...)

Ps = map(Xs) do Xi
    X = Xi[3]
    xs = log10.(getindex.(X, 2))
    ys = getindex.(X, 3)
    zs = getindex.(X, 1)
    kd = kde((xs, ys), weights=zs)
    plot(kd, fill=true, xlim=(-3,-1), ylim=(0,1), colorbar=false)
end 
plot(Ps..., size=(850,650))



## Collapse of local adaptation
# We can try to graph the genome-wide critical m for collapse of local
# adaptation as a function of summary statistics of the barrier?

# Example
k  = 5
Ns̄ = 10
Ls̄ = 1.2
L  = 50
s̄  = Ls̄/L
u  = s̄*0.005
Ne = _Ne2N(Ns̄/s̄, k)
A  = Architecture([randlocus(Exponential(s̄), Beta()) for i=1:L])

Lmx = 1-exp(-Ls̄)
mss = 0:0.05:2.0

ll = map(mss) do ms
    M  = MainlandIslandModel(HapDipDeme(N=N, k=k, A=A, u=u), ms*s̄, ones(L))
    P,_= fixedpointit(M, ones(L))
    p = P[end,:,1]
    pq = P[end,:,2]
    MultilocusIsland.eqload(M, p, pq)
end

plot(ll)
hline!([Lmx])


# Now we try to do something more systematic: do many simulated architectures
# and determine the `m` for which half of the maximal load is obtained.
map(0.01:0.01:0.1) do s̄
    α = 1
    Lmx = 1-exp(-Ls̄)
    A = Architecture([randlocus(Exponential(s̄), Beta(α,α)) for i=1:L])
    function fun(ms) 
        M = MainlandIslandModel(HapDipDeme(N=N, k=k, A=A, u=u), ms*s̄, ones(L))
        P,_= fixedpointit(M, ones(L))
        p = P[end,:,1]
        pq = P[end,:,2]
        MultilocusIsland.eqload(M, p, pq) - Lmx/2
    end
    bisection(fun, 0., 3., 0.02)
end

function bisection(f, a, b, tol)
    fa = f(a)
    fb = f(b)
    c = (a + b) / 2
    y = f(c)
    abs(y) < tol && return c
    (fa < y < fb) && (y < 0) && (return bisection(f, c, b, tol))
    (fa < y < fb) && (y > 0) && (return bisection(f, a, c, tol))
    (fa > y > fb) && (y > 0) && (return bisection(f, c, b, tol))
    (fa > y > fb) && (y < 0) && (return bisection(f, a, c, tol))
end




# s-h correlation
# If we assume that the loci contributing to local adaptation were ate
# mutation-stabilizing selection-drift balance in the ancestral mainland, then
# we might expect higher frequencies for recessive alleles (both phenotypically
# and selectively) at standing variation (s and h will be negatively correlated
# on the mainland). If these are indeed phenotypically recessive, and respond
# to selection on the island, this means that more recessively acting local
# beneficial alleles may contribute, and s and h (of the invading allele) may
# be positively correlated on the island.  We could investigate an ad hoc s by
# h correlated model.

# use two reference points to get a logistic relationship among s and h
function getab(s1, h1, s2, h2)
    b = (log(h2/(1-h2)) - log(h1/(1-h1)))/(log(s2) - log(s1))
    a = log(h1/(1-h1)) - b*log(s1)
    a, b
end

function hfroms(s, a, b, σ)
    lh = a + b*log(s) + rand(Normal(0,σ))
    h  = 1/(1 + exp(-lh))
end

function randlocus2(sd, a, b, σ)
    s = rand(sd)
    h = hfroms(s, a, b, σ)
    DipLocus(-s*h, -s)
end

# p(h) = 1/(h(1-h))*p(lh) = 1/(h(1-h))∫p(lh|s)p(s)ds = 1/(h(1-h))∫N(lh|a + b*log(s),σ)p(s)ds
function marginalh(h, a, b, σ, sd)
    l = log(h/(1-h))
    I, _ = quadgk(s->pdf(sd, s)*pdf(Normal(a + b*log(s), σ), l), 0, 1.)
    return I/(h*(1-h))
end

a, b = getab(0.02, 0.5, 0.1, 0.95)
σ = 1
plot(0:0.001:0.15, map(s->hfroms(s, a, b, σ), 0:0.001:0.15))
plot!(0:0.001:0.15, map(s->hfroms(s, a, b, 0), 0:0.001:0.15))

plot(map(h->marginalh(h, a, b, σ, Exponential(0.02)), 0:0.01:1))

# DFE
s̄  = 0.02
nrep = 100
Xs = tmap([0.1, 0.3, 0.5, 0.7]) do m
    Ps = map([0.5,1.0,1.5]) do Ls
        L = ceil(Int, Ls/s̄)
        sd = Exponential(s̄)
        k  = 5
        Ns = 50.
        N  = _Ne2N(Ns/s̄, k) 
        u  = 0.005*s̄
        plot()
        sss = Float64[]
        hss = Float64[]
        wss = Float64[]
        for i=1:100
            A  = Architecture([randlocus2(sd, a, b, σ) for i=1:L])
            sa = [l.s01 for l in A.loci]
            sb = [l.s11 - 2l.s01 for l in A.loci]
            ss = [l.s11 for l in A.loci]
            hs = [l.s01/l.s11 for l in A.loci]
            M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m*s̄, ones(L))
            P,_= fixedpointit(M, ones(L));
            pm = P[end,:,1];
            ps = singlelocuseq(M);
            ds = abs.(ps .- pm)
            zs = get.(Ref(cs), hs)
            push!(sss, ss...)
            push!(hss, hs...)
            push!(wss, pm...)
        end
        @show m, Ls
        (m, Ls, sss, hss, wss)
    end 
end 


PPs = map(Xs) do X
    map(X) do x
        m, Ls, sss, hss, wss = x
        L = ceil(Int, Ls/s̄)
        binsize = s̄/5
        D = sum(wss) / (L*nrep)
        P1 = stephist(-sss, bins=0:binsize:7s̄, weights=Weights(wss), fill=true,
                      normalize=true, color=:lightgray, xlabel="\$s\$")
        plot!(0:0.001:8s̄, map(x->pdf(Exponential(s̄), x), 0:0.001:8s̄), color=:black, 
              titlefont=9,
              xtickfont=6, ytickfont=6,
              title="\$L\\bar{s}=$Ls, m/\\bar{s} = $m, \\mathbb{E}[D]=$(@sprintf "%.2f" D)\$")
        P2 = stephist(hss, bins=0:0.05:1, weights=Weights(wss), fill=true,
                      normalize=true, color=:lightgray, xlabel="\$h\$")
        plot!(0:0.01:1, map(h->marginalh(h, a, b, σ, Exponential(s̄)), 0:0.01:1), color=:black, xtickfont=6, ytickfont=6)
        plot(P1, P2)
    end
end
PP = deepcopy(vcat(PPs...))
Px = plot(PP..., legend=false, size=(900,500), alpha=0.2, layout=(4,3))

PPs = map(Xs) do X
    map(X) do x
        m, Ls, sss, hss, wss = x
        L = ceil(Int, Ls/s̄)
        xmx = round(quantile(-sss, 0.99), digits=2)
        xbn = range(0, xmx, 20)
        ttl = "\$L\\bar{s} = $Ls, m/\\bar{s} = $m\$"
        P1 = scatter(hss, wss, ms=1, color=:gray, alpha=0.2, 
                ylim=(0,1), xlim=(0,1), title=ttl,
                xlabel="\$h\$", ylabel="\$\\mathbb{E}[p]\$")
        P2 = histogram2d(hss, wss, bins=(0:0.05:1, 0:0.05:1), colorbar=false,
                    xlabel="\$h\$", ylabel="\$\\mathbb{E}[p]\$")
        plot(P1, P2)
    end
end
Px = plot(vcat(PPs...)..., legend=false, size=(1000,700), layout=(4,3), margin=1Plots.mm)

PPs = map(Xs) do X
    map(X) do x
        m, Ls, sss, hss, wss = x
        L = ceil(Int, Ls/s̄)
        xmx = round(quantile(-sss, 0.99), digits=2)
        xbn = range(0, xmx, 20)
        ttl = "\$L\\bar{s} = $Ls, m/\\bar{s} = $m\$"
        zs = get.(Ref(cs), hss)
        P1 = scatter(-sss, wss, ms=1, color=zs, alpha=0.2, 
                title=ttl,
                xlabel="\$s\$", ylabel="\$\\mathbb{E}[p]\$")
    end
end
Px = plot(vcat(PPs...)..., legend=false, size=(700,700), xscale=:log10,
          layout=(4,3), margin=1Plots.mm, xlim=(5e-5, 1e-1))

PPs = map(Xs) do X
    map(X) do x
        m, Ls, sss, hss, wss = x
        L = ceil(Int, Ls/s̄)
        xmx = round(quantile(-sss, 0.99), digits=2)
        xbn = range(0, xmx, 20)
        ttl = "\$L\\bar{s} = $Ls, m/\\bar{s} = $m\$"
        zs = get.(Ref(cs), wss)
        P1 = scatter(-sss, hss, ms=1, color=zs, alpha=0.5, 
                title=ttl,
                xlabel="\$s\$", ylabel="\$h\$")
    end
end
Px = plot(vcat(PPs...)..., legend=false, size=(500,700), xscale=:log10,
          layout=(4,3), margin=1Plots.mm, xlim=(5e-4, 1e-1))


