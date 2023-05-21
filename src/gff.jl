# Gene flow factor calculations
# For the weak migration RV-based argument, the gff is the same for different
# migration models, depending just on the differentiation between the migrant
# pool and the local focal population for which we are computing the gff.

# We have two sorts of cases in terms of numerics/bookkeeping: (1) the case
# where we take an architecture as is and (2) the case with many identical
# unlinked loci, where we compute a gff for each class of loci.

"""
    gff(A, p, pq, y, j)

- `p`  is the allele frequency in the local population
- `pq` is the heterozygosity
- `y`  is the allele frequency in the migrant pool
- `j`  is the locus index for which we are calculating the gff

Note: we may need heterozygosity in the migrant pool as well when we consider
migration at the diploid level?
"""
function gff(A::Architecture, p, pq, y, j)
    @unpack loci, R = A
    x = 0.
    for i=1:length(p)
        i == j && continue
        @unpack s1, s01, s11 = loci[i]
        sa = s1 + s01
        sb = s11 - 2s01
        qi = 1 - p[i]
        x += (sa*(y[i] - qi) + sb*(pq[i] - (1-y[i])*qi)) / R[j,i]
    end
    return exp(x)  # XXX factor 2 disappears here with linkage -> in the rᵢ
end

gff(A::Architecture, p, pq, y) = [gff(A, p, pq, y, j) for j=1:length(A)]


