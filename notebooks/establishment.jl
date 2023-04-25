using ThreadTools, Printf, Plots, PlotThemes, Serialization, Parameters, Random, StatsBase, QuadGK
using MultilocusIsland
theme(:hokusai)

# Simulate until equilibrium
Ls = 1.2
L  = 40
s  = Ls/L
h  = 0.8
N  = 220
k  = 5
A  = Architecture(DipLocus(-s*h, -s), L)

res = map(0.05:0.2:0.85) do ms
    M  = HapDipMainlandIsland(N=N, k=k, m=ms*s, u=0.005*s, arch=A)
    pop, P = simulate(M, 100000)

    # Add a beneficial locus to the architecture
    nrep = 1000
    res  = tmap(1:nrep) do k
        @info k
        eqpop = deepcopy(pop)
        # simulate 1000 generations
        eqpop, P = simulate(M, 1000, pinit=eqpop)
        # modify the architecture
        MM = deepcopy(M)
        push!(MM.arch, DipLocus(s*h, s), 0.5)
        push!(MM.y, 0.)
        push!(eqpop.mainland.genotype, 0)
        # introduce beneficial allele in a single individual
        for (i,x) in enumerate(eqpop.island)
            push!(x, i == N ? 1 : 0)
        end
        # simulate until loss/establishment criterion
        est = false
        n = 0
        while true
            n += 1
            eqpop = MultilocusIsland.generation(MM, eqpop)
            nben = sum(last.(eqpop.island))
            if nben == 0
                break
            elseif nben > N/2
                est = true
                break
            end
        end
        est, n
    end

    estp = sum(first.(res))/nrep

    # calculate me
    E, θ = fixedpointit(M, [1.])
    me = θ[1].m
    gff = me/M.m

    # calculate fixation probability
    function fixp(sa, sb, p0, Ne)
        G(x) = exp(-Ne*x*(2sa + sb*x))
        A, _ = quadgk(G, 0., p0)
        B, _ = quadgk(G, 0., 1.)
        return A/B
    end

    u0 = fixp(s*h - me, s - 2s*h, 1/N, M.Ne)
    u1 = 2s*h - 2me
    ms, estp, u0, u1, me, gff
end

plot( first.(res), getindex.(res, 2), marker=true, label="simulation", xlim=(0,0.9))
plot!(first.(res), getindex.(res, 3), marker=true, label="diffusion")
plot!(first.(res), max.(0., getindex.(res, 4)), marker=true, label="BP", 
      title="\$N_es=6, Ls=$Ls\$", ylabel="\$u\$", xlabel="\$m/s\$", 
      size=(300,250), legend=:topright)


