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
using LogExpFunctions: logit, logistic  
import Random: default_rng
import StatsBase: sample
# these logit/logistic behave better for extreme arguments

include("model.jl")
export HapDipLocus, HapDipMainlandIsland, HapMainlandIsland, sfs, harmonicmean

include("individual-based.jl")
export simulate

include("sampler.jl")
export BetaProposal, BetaSwitchProposal, BetaFlipProposal, gibbs

include("expectation-nlsolve.jl")
export expectedq, expectedsfs

end # module MultilocusIsland
