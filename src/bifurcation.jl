# Deterministic bifurcation analysis for the multilocus model
function tosolve(θ)
	@unpack sa, sb, m, L = θ
	g(p) = exp(2L*p*(sa + sb*(1-p)))
	f(p) = -m*g(p)*p - p*(1-p)*(sa + sb*(1-p))
end

function roots_stab(sa, sb, m, L)
	f  = tosolve((sa=sa, sb=sb, m=m, L=L))
	zs = find_zeros(f, 0.,  1.)
	filter!(!iszero, zs)
	ds = ForwardDiff.derivative.(f, zs)
	zs, ds
end

"""
    findroots_ms(sa, sb, L; stepsize, tol, maxit)

Find the stable and unstable (if any) internal fixed points for a deterministic
multilocus model with `L` equal effect loci. This finds the roots at intervals,
decreasing the step adaptively until they collide.
"""
function findroots_ms(sa, sb, L; stepsize=0.005, tol=1e-2, maxit=10000)
    s = -(sb + 2sa)  # this is the s scale
	ms = 0.
	zs, ds = roots_stab(sa, sb, ms*s, L)
	sol = [(ms, z, d) for (z, d) in zip(zs,ds)]
	i = 0
	while i < maxit && length(zs) >= 1 
		i += 1
		(length(zs) == 2 && abs(zs[1]-zs[2]) < tol) && break
		ms += stepsize
		zs_, ds = roots_stab(sa, sb, ms*s, L)
		if length(zs) == 2 && length(zs_) < 2 
			ms -= stepsize
			stepsize /= 3
		elseif length(zs) == 1 && length(zs_) == 1
			push!(sol, (ms, zs_[1], ds[1]))
            zs = zs_
		elseif length(zs_) == 2
			push!(sol, (ms, zs_[1], ds[1]))
			push!(sol, (ms, zs_[2], ds[2]))
            zs = zs_
		end
	end
	x = first.(sol)
	y = getindex.(sol, 2)
	z = last.(sol)
	stab = z .< 0
	unstab = z .> 0
	x[stab], y[stab], x[unstab], y[unstab] 
end

# The above assumes equal effect loci. We should however be able to generalize
# it to arbitrary effect loci.
