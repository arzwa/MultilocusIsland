"""
    MultilocusIsland

Implements numerical tools for approximating multilocus allele frequency
distributions at equilibrium in the mainland-island model for haploid, diploid
and haplodiplontic populations.
"""
module MultilocusIsland

using StatsBase
using Random
using Printf
using Parameters
using Distributions
using ThreadTools
using QuadGK
using NonlinearSolve
using StaticArrays
using Bijectors
using Bijectors: transform
using LogExpFunctions: logit, logistic  
using ProgressMeter
using Roots, ForwardDiff # for deterministic bifurcation analysis
import Random: default_rng
import StatsBase: sample
import NonlinearSolve: solve
# these logit/logistic behave better for extreme arguments

include("model.jl")
export HapDipMainlandIsland, HapMainlandIsland, Architecture, sfs, harmonicmean
export HapDipLocus, HapLocus, DipLocus

include("ibm.jl")
include("finisland.jl")
export simulate

include("sampler.jl")
include("proposal.jl")
export BetaProposal, BetaSwitchProposal, BetaFlipProposal, gibbs, GibbsSampler, UnitIntervalProposal

include("impliciteq.jl")
export expectedq, expectedsfs, fixedpointit

include("bifurcation.jl")
export findroots_ms, solve

end # module MultilocusIsland
