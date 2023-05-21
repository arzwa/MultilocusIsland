# The exponent of p and (1-p) resp. in Wright's distribution, +1
Afun(N, u, m, pm) = 2N*(u + m*pm)
Bfun(N, u, m, pm) = 2N*(u + m*(1-pm))

# The selection part of Wright's distribution
Cfun(N, p, sa, sb) = exp(-N*p*(2sa + sb*(2-p)))

# The single-locus (log)probability density function
ϕfun(p, N, u, m, pm, sa, sb) = 
    p^(Afun(N,u,m,pm)-1) * (1-p)^(Bfun(N,u,m,pm)-1) * Cfun(N,p,sa,sb)
lϕfun(p, N, u, m, pm, sa, sb) = 
    (Afun(N,u,m,pm)-1)*log(p) + (Bfun(N,u,m,pm)-1)*log(1-p) - N*p*(2sa+sb*(2-p))

function singlelocus_sfs(N, u, m, pm, sa, sb; Δ=0.02, kwargs...)
    Z = log(Zfun((N=N, u=u, m=m, pm=pm, sa=sa, sb=sb); kwargs...))
    ps = (Δ/2):Δ:(1-Δ/2)
    xs = map(p->exp(lϕfun(p, N, u, m, pm, sa, sb) - Z), ps)
    ps, reverse(xs ./ sum(xs))
end

function fpfun(p, θ)
	@unpack sa, sb, N, m, u, pm = θ  # these are the parameters we need
    B = Bfun(N, u, m, pm)
    C = Cfun(N, p, sa, sb)
    -(1-p)^(B-1)*(N * (-p*sb + 2sa + sb*(2-p)) + (B-1)/(1-p)) * C 
end

function gpfun(p, θ)
	@unpack sa, sb, N, m, u, pm = θ
    A = Afun(N, u, m, pm)
    C = Cfun(N, p, sa, sb)
    p^(A-1)*(N*(p*sb - 2sa + sb*(p-2)) + (A-1)/p) * C 
end

function gspfun(p, θ)
	@unpack sa, sb, N, m, u, pm = θ
    A = Afun(N, u, m, pm)
    C = Cfun(N, p, sa, sb)
    p^A*(N*(p*sb - 2sa + sb*(p-2)) + A/p) * C
end

# Normalizing constant of Wright's distribution, using IBP
# could just use Z0fun/Z1fun/Zxfun below (but implemented this one first...)
function Zfun(θ; kwargs...)
	@unpack sa, sb, N, m, u, pm = θ
    A = Afun(N, u, m, pm)
    B = Bfun(N, u, m, pm)
	c = 0.5^(A + B - 1) * exp(-N*(sa + 3sb/4)) * (1/A + 1/B)
	d, _ = quadgk(p->((    p^A)/A)*fpfun(p, θ), 0, 0.5; kwargs...)
	e, _ = quadgk(p->(((1-p)^B)/B)*gpfun(p, θ), 0.5, 1; kwargs...)
	c - d + e
end

# 0 to x integral
function Z0fun(θ, x; kwargs...)
	@unpack sa, sb, N, m, u, pm = θ
    A = Afun(N, u, m, pm)
    B = Bfun(N, u, m, pm)
    a = (1/A)*(x)^A*(1-x)^(B-1)*Cfun(N, x, sa, sb)
	b, _ = quadgk(p->((p^A)/A)*fpfun(p, θ), 0, x; kwargs...)
    return a - b
end

# x to 1 integral
function Z1fun(θ, x; kwargs...)
	@unpack sa, sb, N, m, u, pm = θ
    A = Afun(N, u, m, pm)
    B = Bfun(N, u, m, pm)
    a = (1/B)*(1-x)^B*x^(A-1)*Cfun(N, x, sa, sb)
    b, _ = quadgk(p->(((1-p)^B)/B)*gpfun(p, θ), x, 1; kwargs...)
    return a - b
end

# General x to y integral
function Zxfun(θ, x, y; kwargs...)
	@unpack sa, sb, N, m, u, pm = θ
    x ≈ 0 && y ≈ 1 && return Zfun(θ; kwargs...)
    x ≈ 0 && return Z0fun(θ, y; kwargs...)
    y ≈ 1 && return Z1fun(θ, x; kwargs...)
    I, _ = quadgk(p->ϕfun(p, N, u, m, pm, sa, sb), x, y; kwargs...)
    return I
end

function Yfun(θ; kwargs...)
	@unpack sa, sb, N, m, u, pm = θ
    A = Afun(N, u, m, pm)
    B = Bfun(N, u, m, pm)
	C(p) = Cfun(N, p, sa, sb)
	L, _ = quadgk(p->p^A * (1-p)^(B-1) * C(p), 0., 0.5; kwargs...)
	R1 = (1/B)*(1/2)^(B+A)*C(0.5)
	R2, _ = quadgk(p->(((1-p)^B)/B)*gspfun(p, θ), 0.5, 1; kwargs...)
	return L + R1 + R2
end

function Epqfun(θ; kwargs...)
	@unpack sa, sb, N, m, u, pm = θ
    A = Afun(N,u,m,pm)
    B = Bfun(N,u,m,pm)
    N, _ = quadgk(p -> p^A * (1-p)^B * Cfun(N,p,sa,sb), 0., 1.)
    Z = Zfun(θ; kwargs...)
    return N/Z
end

# Expected value
Epfun(θ; kwargs...) = Yfun(θ; kwargs...) / Zfun(θ; kwargs...)

