using MultilocusIsland, StatsBase, ColorSchemes, ThreadTools, Printf, QuadGK
using Plots, PlotThemes, Distributions, Parameters; theme(:hokusai)
using LogExpFunctions
cs = ColorSchemes.viridis

function singlelocuseq(M)
    map(1:length(M.D.A)) do i
        DD = reconstruct(M.D, A=M.D.A[i:i])
        MM = MainlandIslandModel(DD, M.m1)
        P,_= fixedpointit(MM, ones(1))
        P[end,:,1]
    end |> x->first.(x)
end

# amounts more or less to the same I suppose.
function msedfe1(dfe, L, N, k, u, m; n=100) 
    sss = Float64[]
    hss = Float64[]
    wss = Float64[]
    dss = Float64[]
    for i=1:n
        A  = Architecture([randlocus(dfe) for i=1:L])
        sa = [l.s01 for l in A.loci]
        sb = [l.s11 - 2l.s01 for l in A.loci]
        ss = [l.s11 for l in A.loci]
        hs = [l.s01/l.s11 for l in A.loci]
        M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m)
        P,_= fixedpointit(M, ones(L));
        pm = P[end,:,1];
        ps = singlelocuseq(M);
        ds = abs.(ps .- pm)
        push!(sss, ss...)
        push!(hss, hs...)
        push!(wss, pm...)
        push!(dss, ds...)
    end
    sss, hss, wss, dss
end 

Ns  = 20
s   = 0.01
Ne  = Ns/s
k   = 5
N   = _Ne2N(Ne, k)
u   = s * 0.005
κs  = [4, 1, 1/4]
λ   = 1/s
Ls  = 1.0
L   = round(Int, Ls/s)
mms = 0.05:0.15:0.65

res = map(κs) do κ
    tmap(mms) do ms
        @show ms, κ
        dfe1 = IndependentDFE(Gamma(κ, s/κ), Beta(1,1))
        msedfe1(dfe1, L, N, k, u, ms*s)
    end
end

function scat(out, κs, mms; yh=(0,Inf), xmx=0.052)
    Ps = map(zip(κs, out)) do (κ, x)
        map(zip(mms, x)) do (m, y)
            ss, hs, ws, ds = y    
            yl = m == 0.05
            xl = κ == 0.25
            idx = sample(1:length(ws), 1000, replace=false)
            ss = log10.(-ss)
            ss_ = (ss .- minimum(ss)) ./ (maximum(ss) - minimum(ss))
            P1 = scatter(ws[idx] .- ds[idx], ws[idx], 
                    color=get.(Ref(cs), hs[idx]),
                    alpha=0.5, ms=3,
                    title="\$\\kappa = $κ, m/\\bar{s} = $m\$", 
                    legend=false, xlim=(0,1), ylim=(0,1),
                    ylabel= yl ? "\$\\mathbb{E}[p_i]\$ multilocus" : "",
                    xlabel= xl ? "\$\\mathbb{E}[p_i]\$ single locus" : "")
            plot!(P1, x->x, color=:gray)
            if κ==0.25 &&  m==0.65
                plot!(P1, repeat([0.9], 20), range(0.1,0.8,20), 
                      lw=4, color=get.(Ref(cs), range(0,1,20)))
                annotate!(0.82, 0.45, text("\$h\$", 9, :left))
                annotate!(0.93, 0.1, text("\$0\$", 9, :left))
                annotate!(0.93, 0.8, text("\$1\$", 9, :left))
            end
#            lws = log10.(ws)
#            xss = [0.01, 0.1, 1]
            P2 = stephist(ws, fill=true, color=:black, alpha=0.2,
                          legend=false, norm=true, framestyle=:default, 
                          yshowaxis=false, xlim=(0,1))
                          #xlim=(-2.2,0),
                          #xticks=(log10.(xss), string.(xss)))
            plot(P1, P2, layout=grid(2,1, heights=[0.7,0.3]))
        end
    end
    Ps = vcat(Ps...)
#    annotate!(Ps[1], -0.35, 1.105, text("\$\\mathrm{(A)}\$", 10))
    # colorbar
    Ps = plot(Ps..., size=(890,850), margin=2Plots.mm, titlefont=9,
              layout=(3,5))
end
scat(res, κs, mms)

function dens(out, κs, mms; yh=(0,Inf), xmx=0.052)
    Ps = map(zip(κs, out)) do (κ, x)
        map(zip(mms, x)) do (m, y)
            ss, hs, ws, ds = y    
            yl = m == 0.05
            xl = κ == 0.25
            DD = @sprintf "\$\\bar{\\Delta} = %.2f\$" sum(ws) / length(ws)
            idx = sample(1:length(ws), 1000, replace=false)
            stephist(ws, bins=0:0.025:1, fill=true, norm=:probability,
                     color=:black, alpha=0.2, linealpha=0, ylim=(0,0.33),
                     title="\$\\kappa = $κ, m/\\bar{s} = $m\$", 
                    legend=false, xlim=(0,1),
                    xlabel= xl ? "\$\\mathbb{E}[p_i]\$ multilocus" : "",
                    ylabel= yl ? "frequency" : "")
            annotate!(0.1, 0.3, text(DD, 8, :left))
#            stephist!([w for (w,h) in zip(ws, hs) if h < 1/3], color=1,
#                      bins=0:0.025:1, norm=:probability, fill=true, alpha=0.5)
#            stephist!([w for (w,h) in zip(ws, hs) if h > 2/3], color=3,
#                      bins=0:0.025:1, norm=:probability, fill=true, alpha=0.5)
        end
    end
    Ps = vcat(Ps...)
    Ps = plot(Ps..., size=(890,520), margin=2Plots.mm, titlefont=9,
              layout=(3,5))
end
dens(res, κs, mms)
