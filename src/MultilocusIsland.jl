"""
    MultilocusIsland

Implements numerical tools for approximating multilocus allele frequency
distributions at equilibrium in the mainland-island model for haploid, diploid
and haplodiplontic populations.
"""
module MultilocusIsland

using StatsBase
using Random
using Parameters
using Distributions
using ThreadTools
using QuadGK
using NonlinearSolve
using StaticArrays
using Bijectors
using LogExpFunctions: logit, logistic  
import Random: default_rng
import StatsBase: sample
# these logit/logistic behave better for extreme arguments

include("model.jl")
export HapDipLocus, HapDipMainlandIsland, HapMainlandIsland, Architecture, sfs, harmonicmean

include("ibm.jl")
export simulate

include("sampler.jl")
include("proposal.jl")
export BetaProposal, BetaSwitchProposal, BetaFlipProposal, gibbs, GibbsSampler, UnitIntervalProposal

include("impliciteq.jl")
export expectedq, expectedsfs

end # module MultilocusIsland
