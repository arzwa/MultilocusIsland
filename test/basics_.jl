using MultilocusIsland

L = 50
Ls = 1.
s = Ls/L
Ns = 8
k = 5
N = _Ne2N(Ns/s, k)
u = s*0.005
m = 0.2*s
r = 10s

import MultilocusIsland: MainlandIslandModel, HapDipDeme
A = Architecture(HapLocus(-s), L, r) 
M = MainlandIslandModel(HapDipDeme(N=N, k=k, u=u, A=A), m, ones(L)) 
_, P = simulate(M, 51000, zeros(L), drop=1000, thin=2)
Q, _ = MultilocusIsland.fixedpointit_linkage(M, ones(L));

histogram(vec(P), bins=0:0.025:1)
vline!([1 .- Q[end,:,1]])

# Haplodiplontic model
A = Architecture([HapDipLocus(-s*(1-τ), -s*h*τ, -s*τ) for i=1:L], fill(0.5, L)) 
M = HapDipMainlandIsland(N=N, k=k, m=m, arch=A, u=u)
_, P = simulate(M, 1000)
G = GibbsSampler([UnitIntervalProposal() for i=1:L])
Q, l = gibbs(M, G, rand(L), 10000)

fp = MultilocusIsland._fixedpointit(M, [0.5])
# Note that the Gibbs sampler/fixed point iteration for the haploid model is
# implemented as a special case of the diploid model (this is not true for the
# IBM).

A = Architecture([HapDipLocus(-s*(1-τ), -s*h*τ, -s*τ) for i=1:L], fill(0.5, L)) 
M = HapDipMainlandIsland(N=N, k=k, m=m, arch=A, u=u)
fp = MultilocusIsland.fixedpointit(M, [0.5])

