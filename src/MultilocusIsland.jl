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

include("architecture.jl")
export Architecture, HapDipLocus, HapLocus, DipLocus

include("models.jl")
export simulate, MainlandIslandModel, FiniteIslandModel, HapDipDeme

#include("sampler.jl")
#include("proposal.jl")
#export BetaProposal, BetaSwitchProposal, BetaFlipProposal, gibbs, GibbsSampler, UnitIntervalProposal
#include("gff.jl")
#include("proposal.jl")
#include("_mrfsampler.jl")

include("integration.jl")
include("fixedpointit.jl")
export expectedq, expectedsfs, fixedpointit, summarize_arch

include("bifurcation.jl")
export findroots_ms, solve

include("utils.jl")
export _Ne2N, sfs, harmonicmean

include("dfe.jl")
export IndependentDFE, CKExponential, CKGamma, Logisticsbyh, randlocus

end # module MultilocusIsland
