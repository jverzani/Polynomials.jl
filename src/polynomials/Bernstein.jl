

struct Bernstein{T <: Number} <: Polynomials.AbstractPolynomial{T}
N::Int
coeffs::Vector{T}
var::Symbol # would like to relax this
function Bernstein{T}(N::Int, nu::Int, var::Symbol=:x) where {T <: Number}
    ## Basis vectors b_{\nu,N}
    0 <= nu <= N || throw(ArgumentError("Must have 0 ≤ nu ≤ N"))
    cs = zeros(T, N+1)
    cs[nu+1] = one(T)
    new{T}(N, cs, var)
end
function Bernstein{T}(coeffs::AbstractVector{T}, var::Symbol=:x) where {T <: Number}
    isempty(coeffs)  && return new{T}(0, zeros(T, 0), var)
    last_nz = findlast(!iszero, coeffs)
    last = max(1, last_nz === nothing ? 0 : last_nz)
    N =   last - 1
    return new{T}(N,  coeffs[1:last], var)
end
end

Polynomials.@register Bernstein

function Polynomials.showterm(io::IO, ::Type{Bernstein{T}}, pj::T, var, j, first::Bool, mimetype) where {T}
    iszero(pj) && return true
    !first && print(io, " + ")
    print(io, "$pj ⋅ β(N, $j)($var) ")
    return true
end


function Base.convert(P::Type{<:Polynomial{T}}, ch::Bernstein) where {T}
    N = ch.N
    out = P(zeros(T,1), ch.var)
    x = P([zero(T),one(T)], ch.var)
    @inbounds for (i,ai) in enumerate(coeffs(ch))
        nu = i - 1
        out += ai *  binomial(N, nu) * x^nu * (1-x)^(N-nu)
    end
    out

end

function Base.convert(C::Type{<:Bernstein{T}}, p::Polynomial) where {T}

    N = degree(p)
    cs = zeros(T, N+1)

    for (i, a_nu) in enumerate(coeffs(p))
        k = i - 1
        nk = binomial(N,k)
        for j in k:N
            cs[j+1] += a_nu/nk*binomial(j,k)
        end
    end
    Bernstein{T}(cs, p.var)
end

domain(::Type{<:Bernstein}) = Polynomials.Interval(0, 1)


function (ch::Bernstein{T})(x::S) where {T,S}
    # XXX make more efficient
    x ∉ domain(ch) && error("$x outside of domain")
    R = promote_type(T, S)
    length(ch) == 0 && return zero(R)
    convert(Polynomial{T}, ch)(x)
end


function integrate(p::Bernstein{T}, C::S) where {T,S <: Number}
    R = promote_type(eltype(one(T) / 1), S)
    n = p.N
    N = n + 1
    cs = zeros(R, N+1)

    x = variable(p)
    @inbounds for (i,a_nu) in enumerate(coeffs(p))
        nu = i - 1
        for j in (nu+1):N
            cs[j+1] += a_nu/(n+1)
        end
    end
    Bernstein{R}(cs, p.var)
end




function derivative(p::Bernstein{T}, order::Integer = 1) where {T}
    N =  p.N  - 1
    cs = coeffs(p)
    N < 0 &&  return p
    order ==  0 && return p
    csp = zeros(eltype(cs), N+1)
    # nu = 0
    for (i,a_nu) in enumerate(cs)
        nu = i - 1
        nu > 0  && (csp[nu] += a_nu *  (N+1))
        nu <= N && (csp[nu+1] -= a_nu *  (N+1))
    end
    csp
    pp = Bernstein{T}(csp, p.var)
    derivative(pp, order-1)
end


function Base.:+(p1::Bernstein{T}, p2::Bernstein{S}) where {T, S}

    p1.var == p2.var || throw(ArgumentError("p1 and p2 must have the same symbol"))

    R = promote_type(T, S)
    N1, N2 = p1.N, p2.N

    if N1 == N2
        return Bernstein{R}(coeffs(p1) +coeffs(p2), p1.var)
    else
        q1 = convert(Polynomial{R}, p1)
        q2 = convert(Polynomial{R}, p2)
        return convert(Bernstein{R}, q1+q2)
    end
end

function Base.:*(p1::Bernstein{T}, p2::Bernstein{S}) where {T,S}

    p1.var == p2.var || throw(ArgumentError("p1 and p2 must have the same symbol"))

    R = promote_type(T, S)
    q1 = convert(Polynomial{T}, p1)
    q2 = convert(Polynomial{S}, p2)
    return convert(Bernstein{R}, q1*q2)

end

function Base.divrem(num::Bernstein{T}, den::Bernstein{S}) where {T,S}
    ## XXX
    p1 = convert(Polynomial{T}, num)
    p2 = convert(Polynomial{S}, den)
    q,r = divrem(p1, p2)
    R = eltype(q)

    convert.(Bernstein{R}, (q,r))
end

function companion(p::Bernstein{T}) where {T}
    companion(convert(Polynomial{T}, p))
end
