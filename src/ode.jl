# Deterministic model
using SteadyStateDiffEq, DifferentialEquations

# XXX not updated nor tested
function deterministic_eq(M::MainlandIslandModel, p0=1.)
    θ = summarize_arch(M)
    f(p, θ) = ode!(similar(p), p, θ, NaN)
    p0 = fill(p0, θ.K)
    prob = NonlinearProblem{false}(f, p0, θ)
    return solve(prob, DynamicSS(Rodas5P()))
end

function solve_ode(M::MainlandIslandModel, t, p0=1.)
    θ = summarize_arch(M)
    X = ODEProblem(ode!, fill(p0, θ.K), (0,t), θ)
    return solve(X)
end

# for solving with DifferentialEquations
function ode!(dp, p, θ, t)
	@unpack m, loci, L, γ, y = θ
	xs = map(zip(p, loci, γ, y)) do (pj, lj, wj, yj)
		@unpack s1, s01, s11 = lj
        aj = s1 + s01 + (s11 - 2s01)*(1-pj) 
        aj, wj*(yj-(1-pj))*aj
	end
	a = first.(xs)
	g = exp(2L*sum(last.(xs)))  # gene flow factor
	for i=1:length(p)
        dp[i] = -m*g*(y[i]-(1-p[i])) - p[i]*(1-p[i])*a[i]
	end
    return dp
end

