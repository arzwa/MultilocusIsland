using MultilocusIsland, StatsBase, ColorSchemes
using Plots, PlotThemes, Distributions, Parameters; theme(:hokusai)
cs = ColorSchemes.viridis
Ne2N(Ne, k) = ceil(Int, Ne/(2k) + Ne)

function randlocus(sd, hd)
    s = rand(sd)
    h = rand(hd)
    DipLocus(-s*h, -s)
end

# deterministic single-locus predictions
function singlelocuseq_det(M)
    map(1:length(M.arch)) do i
        @unpack s01, s11 = M.arch[i]
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
    map(1:length(M.arch)) do i
        MM = reconstruct(M, arch=M.arch[i:i], y=M.y[i:i])
        P,_= fixedpointit(MM, ones(1))
        P[end,:,1]
    end |> x->first.(x)
end

function standardize(xs)
    xmn, xmx = extrema(xs)
    xs .-= xmn
    xs ./ (xmx - xmn)
end

L = 100
PPs = map([0.05, 0.1, 0.2]) do m
    Ps = map([0.5,1.0,1.5]) do Ls
        s̄ = Ls/L
        #Ls = 1.0
        #L  = ceil(Int, Ls/s̄)
        sd = Exponential(s̄)
        hd = Beta()
        k  = 5
        Ns = 8.  
        N  = Ne2N(Ns/s̄, k) 
        u  = 0.005*s̄
        plot()
        for i=1:5
            A  = Architecture([randlocus(sd, hd) for i=1:L])
            sa = [l.s01 for l in A.loci]
            sb = [l.s11 - 2l.s01 for l in A.loci]
            hs = [l.s01/l.s11 for l in A.loci]
            M  = HapDipMainlandIsland(N=N, k=k, m=m*s̄, u=u, arch=A)
            xs = summarize_arch(M)
            P,_= fixedpointit(M, ones(xs.K));
            pm = P[end,:,1];
            ps = singlelocuseq(M);
            zs = get.(Ref(cs), hs)
            scatter!(ps, pm, color=zs, markerstrokecolor=zs, ms=3, alpha=0.8)
        end
        plot!(x->x, color=:black, xlim=(0,1), ylim=(0,1),
              titlefont=9,
              xtickfont=7, ytickfont=7,
              title="\$L\\bar{s}=$Ls, m/\\bar{s} = $m\$", 
              xlabel=m == 0.2 ? "singlelocus" : "", 
              ylabel=Ls == 0.5 ? "multilocus" : "")
    end
    vcat(Ps...)
end 

PP = deepcopy(vcat(PPs...))
plot!(PP[end], repeat([0.9], 20), range(0.1,0.8,20), 
      lw=4, color=get.(Ref(cs), range(0,1,20)))
annotate!(0.82, 0.45, text("\$h\$", 9, :left))
annotate!(0.93, 0.1, text("\$0\$", 9, :left))
annotate!(0.93, 0.8, text("\$1\$", 9, :left))
plot(PP..., legend=false, size=(600,650))

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
M  = HapDipMainlandIsland(N=N, k=k, m=m*s̄, u=u, arch=A)
xs = summarize_arch(M)
P,_= fixedpointit(M, ones(xs.K));
p  = P[end,:,1]
pq = P[end,:,2]

ys = expectedsfs(M, p, pq, step=0.005, f=log10)
plot(ys)

_,Q = simulate(M, 210000, drop=10000, thin=2)

PP = plot(title="\$L=$L, Ls = $Ls, N_e s = $Ns, m/s = $m\$", size=(450,320))
map(enumerate(1:5)) do (c,i)
    h = A[i].s01/A[i].s11
    s = -A[i].s11
    plot!(ys[i], color=c, label=@sprintf "locus %d \$(s=%.3f, h=%.2f)\$" i s h)
    scatter!(sfs(Q[:,i], step=0.05, f=log10), markerstrokecolor=c, color=c, label="")
end
plot(PP, legend=:bottomleft, ylabel="\$\\log_{10}\\phi(p)\$", xlabel="\$p\$")


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
