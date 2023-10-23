
_Ne2N(Ne, k) = ceil(Int, Ne/(2k) + Ne)

harmonicmean(x, y) = 1/(1/x + 1/y)

function lognormalize(l)
    x = exp.(l .- maximum(l))
    return x ./ sum(x)
end

function sfs(ps; step=0.02, f=identity, density=true)
    hs = fit(Histogram, ps, 0:step:1+step) 
    hs = StatsBase.normalize(hs, mode=:probability)
    es = collect(hs.edges)[1][2:end-1] .- step/2
    ws = hs.weights[1:end-1]    # make sure we add the fixed states to the last bin 
    ws[end] += hs.weights[end] 
    if density
        ws ./= step
    end
    es, f.(ws)
end

# Haldane's mapping function (map distance -> P recombination)
haldane(y) = 0.5*(1-exp(-y))

# the inverse of Haldane's mapping function (P rec -> map distance)
invhaldane(x) = -log(1 - 2x)

# Compute recombination rates across a linear genome based on pairwise
# recombination fractions
function rrates(xs, j)
    ys = invhaldane.(xs)
    left = reverse(cumsum(ys[j-1:-1:1]))
    rght = cumsum(ys[j:end-1])
    [haldane.(left); NaN; haldane.(rght)]
end

# recombination rate matrix
rrates(xs) = hcat(rrates.(Ref(xs), 1:length(xs))...)

