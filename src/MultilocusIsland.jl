"""
    MultilocusIsland
"""
module MultilocusIsland

using Printf
using QuadGK
using Random
using Roots
using ForwardDiff
using Reexport
using Parameters
using StatsBase
using Distributions
using LinearAlgebra
using ProgressMeter
@reexport using WrightDistribution
import WrightDistribution: expectedpq

include("architecture.jl")
include("recombination.jl")
include("deme.jl")
include("mainlandisland.jl")
include("metapopulation.jl")
include("simulation.jl")
include("effectivemigration.jl")
include("utils.jl")
include("deterministic.jl")

vvcat(x) = vcat(x...)

export Locus, Architecture, Deme, FiniteIslands, MainlandIsland, FixedMainland
export generation!, initpop, vvcat, sasb, sehe, Ne2N, eqpdf, simulate!
export humanmap, humanmap2, flymap, flymap2

end
