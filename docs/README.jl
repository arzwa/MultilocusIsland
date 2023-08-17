# # MultilocusIsland.jl

# This julia package implements:
#
# 1. Individual-based simulations of a multilocus, biallelic, mainland-island
#    population genetics model for a haplodiplontic/haploid/diploid life cycle.
# 2. Numerical tools to calculate effective migration rates and approximate
#    allele frequency distributions in the above model using diffusion theory.
#
# To install the package:
#
# 1. Install `julia` from [https://julialang.org/](https://julialang.org/).
# 2. Open up a terminal and start a REPL by entering `julia` in the prompt.
# 3. Type `]` to enter the julia package manager.
# 4. Type `add https://github.com/arzwa/MultilocusIsland` to install the package.
#
# After executing these steps, one should be able to load the package using
# `using MultilocusIsland` in the julia REPL.

# ## Examples

using MultilocusIsland, Plots


