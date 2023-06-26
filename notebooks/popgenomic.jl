# Simulate population genomic data. Eventually we'd simulate doing full forward
# in time simulation of an evolutionary scenario I guess (maybe using Slim),
# but at first it would be good to just simulate according to thegenerative
# model and do inference. 
#
# What's the generative model? 
# - We assume a WF population
# - We have some genetic architecture of local adaptation
# - We have a migration rate m
# - We assume migration-selection-mutation-drift balance for the selected
#   sites, yielding L MSPs
# - We assume neutral sites linked to the MSPs
# - We generate samples for the neutral sites according to a structured
#   coalescent
# 
# Well, we can just simulate this forward in time using our current tools. So
# we can directly focus on the inferential question, which involves of course
# simulation from the above generative model in an ABC framework...
# So we can simulate this model forward in time, the question is, to what
# extent can we approximate it to the me-structured coalescent thing and to
# what extent can we estimate parameters...

# The only piece that is really missing for this is structured coalescent
# simulations...
function randlocus(sd, hd, p0)
    rand() < p0 && return DipLocus(0., 0.)
    s = rand(sd)
    h = rand(hd)
    DipLocus(-s*h, -s)
end

k  = 5
Ns̄ = 10
Ls̄ = 1.2
p0 = 0.9
L  = 1000
s̄  = Ls̄/(L*(1-p0))
u  = s̄*0.005
ms = 0.2
Ne = _Ne2N(Ns̄/s̄, k)
A  = Architecture([randlocus(Exponential(s̄), Beta(), p0) for i=1:L], (rand(L) .^2) ./ 2)
y  = [l.s11== 0. ? 0.5 : 1. for l in A.loci]
M  = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), ms*s, y)

P, _ = MultilocusIsland.fixedpointit_linkage(M, ones(L));
p  = P[end,:,1]
pq = P[end,:,2]
mes = MultilocusIsland.gffall(M, p, pq)

plot(getindex.(mes, Ref(:m)))

