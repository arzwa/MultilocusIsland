using Random, MultilocusIsland
using Plots, PlotThemes; theme(:hokusai)

N  = 200
k  = 5
Ne = harmonicmean(N, 2N*k)
L  = 40 
s  = 4/Ne
h  = 0.8
τ  = 0.1

mss = 0:0.1:1
rs = [0.01, 0.05, 0.1, 0.25, 0.5]
results = map(rs) do r
    p = map(mss) do ms
        @info ms, r
        A = Architecture([HapDipLocus(-s*(1-τ), -s*h*τ, -s*τ) for i=1:L], fill(r, L)) 
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, arch=A, u=s*0.005)
        _, P = simulate(M, 11000)
        mean(P[1001:end,:])
    end
    (N=N, k=k, L=L, h=h, τ=τ, r=r, s=s, m=mss, p=p)
end

ker = BetaFlipProposal(0.5, 1.0, 0.1)
smples = map(results) do x
    p = map(x.m) do ms
        @info (x, ms)
        A = Architecture([HapDipLocus(-s*(1-τ), -s*h*τ, -s*τ) for i=1:L], fill(x.r, L)) 
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, arch=A, u=s*0.005)
        Q, _ = gibbs(M, ker, rand(L), 10000+1000, drop=1000)
        mean(Q)
    end
end

plot()
map(enumerate(results)) do (i,x)
    plot!(x.m, 1 .- x.p, label="\$r=$(x.r)\$", color=i, ms=3,
          markerstrokecolor=i, marker=true)
    pp = smples[i]
    plot!(x.m, 1 .- pp, label="\$r=$(x.r)\$", color=i, ms=3,
          markerstrokecolor=i, marker=true, ls=:dot)
end
plot!()


# an example to focus on
L = 20
s = 0.01
h = 0.8
τ = 0.0
r = s*1
N = 220
k = 5
u = s*0.005

xs = map(1:5) do rep
    A = Architecture([HapDipLocus(-s*(1-τ), -s*h*τ, -s*τ) for i=1:L], fill(r, L)) 
    map(mss) do ms
        @info (rep, ms)
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, arch=A, u=u)
        _, P = simulate(M, 11000, drop=1000, thin=10)
        vec(mean(P, dims=1))
    end
end

#ker = BetaFlipProposal(0.5, 1.0, 0.1)
ker = BetaProposal(0.2)
ys = map(1:5) do rep
    A = Architecture([HapDipLocus(-s*(1-τ), -s*h*τ, -s*τ) for i=1:L], fill(r, L)) 
    map(mss) do ms
        @info (rep, ms)
        M = HapDipMainlandIsland(N=N, k=k, m=ms*s, arch=A, u=s*0.005)
        Q, _ = gibbs(M, ker, rand(L), 5100, drop=100)
        vec(mean(Q, dims=1))
    end
end

plot()
xmean = mean(xs)
ymean = mean(ys)
map(enumerate(1:4:4*3)) do (i,j)
    plot!(mss, 1 .- getindex.(xmean, j), color=i, ms=3, markerstrokecolor=i,
          marker=true)
    plot!(mss, 1 .- getindex.(ymean, j), color=i, ms=3, markerstrokecolor=i,
          marker=true, ls=:dot)
end
plot!()

plot( 1 .- mean.(xmean))
plot!(1 .- mean.(ymean))

