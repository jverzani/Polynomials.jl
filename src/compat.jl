## Leave these deprecations for later
poly(r, var = :x) = fromroots(Polynomial, r; var = var)
polyval(p::AbstractPolynomial, x::Number) = p(x)
polyval(p::AbstractPolynomial, x) = p.(x)

polyint(p::AbstractPolynomial, k = 0) = integral(p, k)
polyint(p::AbstractPolynomial, a, b) = integrate(p, a, b)
polyder(p::AbstractPolynomial, ord = 1) = derivative(p, ord)
polyfit(x, y, n = length(x) - 1) = fit(Polynomial, x, y; deg = n)
polyfit(x, y, sym::Symbol) = fit(Polynomial, x, y; var = sym)
Poly(as, var=:x) = Polynomial(as, var)

padeval(PQ::Pade, x::Number) = PQ(x)
padeval(PQ::Pade, x) = PQ.(x)

## These are to be deprecated
export Poly, poly, polyval, polyint, polyder, polyfit, padeval

#

function Base.getproperty(p::AbstractPolynomial, nm::Symbol)
    if nm == :a
        Base.depwarn(
            "AbstractPolynomial.a is deprecated, use AbstracPolynomial.coeffs instead.",
            Symbol("Base.getproperty"),
        )
        return getfield(p, :coeffs)
    end
    return getfield(p, nm)
end

Base.zero(::Type{Poly{T}}) where {T} = Polynomial(zeros(T,1))
Base.one(::Type{Poly{T}}) where {T} = Polynomial(ones(T,1))
