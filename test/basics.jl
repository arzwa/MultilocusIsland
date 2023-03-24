using MultilocusIsland

L = 20
s = 0.04
h = 0.5
τ = 1.
N = 220
k = 5
u = s*0.005
m = 0.01

# Haploid model
A = Architecture([HapLocus(-s) for i=1:L], fill(0.5, L)) 
M = HapMainlandIsland(N=N, m=m, arch=A, u=u)
_, P = simulate(M, 50000)

# Haplodiplontic model
A = Architecture([HapDipLocus(-s*(1-τ), -s*h*τ, -s*τ) for i=1:L], fill(0.5, L)) 
M = HapDipMainlandIsland(N=N, k=k, m=m, arch=A, u=u)
_, P = simulate(M, 1000)
G = GibbsSampler([UnitIntervalProposal() for i=1:L])
Q, l = gibbs(M, G, rand(L), 10000)

fp = fixedpointit(M, [0.5])
# Note that the Gibbs sampler/fixed point iteration for the haploid model is
# implemented as a special case of the diploid model (this is not true for the
# IBM).

