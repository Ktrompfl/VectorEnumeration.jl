using Nemo
using VectorEnumeration
using BenchmarkTools


#E₁
field = GF(2)
algebra, (s, t, r, s′, t′, r′) = free_associative_algebra(field, ["s", "t", "r", "s'", "t'", "r'"])
T = elem_type(algebra)
inverses = Union{T, Tuple{T, T}}[(s, s′), (t, t′), (r, r′)]
relators = Union{T, Tuple{T, Weight}}[ t′*r*t*r′^2, r′*s*r*s′^2, s′*t*s*t′^2 ]
subgens = NTuple{1, T}[]

println("vectorenumeration:")
@btime vector_enumeration(algebra, inverses, relators, subgens)

idealgen = [s*s′ - 1, t*t′ - 1, r*r′ - 1, t′*r*t*r′^2 - 1, r′*s*r*s′^2 - 1, s′*t*s*t′^2 - 1 ]

println("groebnerbasis:")
@btime Generic.groebner_basis(idealgen)

#Cox
field = GF(2)
algebra, (a, b, a′, b′) = free_associative_algebra(field, ["a", "b", "a'", "b'"])
T = elem_type(algebra)
inverses = Union{T, Tuple{T, T}}[(a, a′), (b, b′)]
relators = Union{T, Tuple{T, Weight}}[ a^6, b^6, (a*b)^2, (a^2*b^2)^2, (a^3*b^3)^5 ]
subgens = NTuple{1, T}[]

println("vectorenumeration:")
@btime vector_enumeration(algebra, inverses, relators, subgens)

idealgen = [a*a′ - 1, b*b′ - 1, a^6 - 1, b^6 - 1, (a*b)^2 - 1, (a^2*b^2)^2 - 1, (a^3*b^3)^5 - 1]

println("groebnerbasis:")
@btime Generic.groebner_basis(idealgen)

#B₂,₄
field = GF(2)
algebra, (a, b, a′, b′) = free_associative_algebra(field, ["a", "b", "a'", "b'"])
T = elem_type(algebra)
inverses = Union{T, Tuple{T, T}}[(a, a′), (b, b′)]
relators = Union{T, Tuple{T, Weight}}[ a^4, b^4, (a*b)^4, (a′*b)^4, (a^2*b)^4, (a*b^2)^4, (a^2*b^2)^4, (a′*b*a*b)^4, (a*b′*a*b)^4 ]
subgens = NTuple{1, T}[]

println("vectorenumeration:")
@btime vector_enumeration(algebra, inverses, relators, subgens)

idealgen = [a*a′ - 1, b*b′ - 1, a^4 - 1, b^4 - 1, (a*b)^4 - 1, (a′*b)^4 - 1, (a^2*b)^4 - 1, (a*b^2)^4 - 1, (a^2*b^2)^4 - 1, (a′*b*a*b)^4 - 1, (a*b′*a*b)^4 - 1]

println("groebnerbasis:")
@btime Generic.groebner_basis(idealgen)

#S₇
field = GF(2)
algebra, (a, b, a′, b′) = free_associative_algebra(field, ["a", "b", "a'", "b'"])
T = elem_type(algebra)
inverses = Union{T, Tuple{T, T}}[(a, a′), (b, b′)]
relators = Union{T, Tuple{T, Weight}}[ a^7, b^2, (a*b)^6, (a*b*a′*b′)^3, (a*a*b*a′*a′*b′)^2, (a*a*a*b*a′*a′*a′*b′)^2 ]
subgens = NTuple{1, T}[]

println("vectorenumeration:")
@btime vector_enumeration(algebra, inverses, relators, subgens)

idealgen = [a*a′ - 1, b*b′ - 1, a^7 - 1, b^2 - 1, (a*b)^6 - 1, (a*b*a′*b′)^3 - 1, (a*a*b*a′*a′*b′)^2 - 1, (a*a*a*b*a′*a′*a′*b′)^2 - 1]

println("groebnerbasis:")
@btime Generic.groebner_basis(idealgen)

#PSL₂(11)
field = GF(2)
algebra, (a, b, a′, b′) = free_associative_algebra(field, ["a", "b", "a'", "b'"])
T = elem_type(algebra)
inverses = Union{T, Tuple{T, T}}[(a, a′), (b, b′)]
relators = Union{T, Tuple{T, Weight}}[ a^11, b^2, (a*b)^3, (a^4*b*a′^5*b)^2 ]
subgens = NTuple{1, T}[]

println("vectorenumeration:")
@btime vector_enumeration(algebra, inverses, relators, subgens)

idealgen = [a*a′ - 1, b*b′ - 1, a^11 - 1, b^2 - 1, (a*b)^3 - 1, (a^4*b*a′^5*b)^2 - 1]

println("groebnerbasis:")
@btime Generic.groebner_basis(idealgen)

#PSL₂(17)
field = GF(2)
algebra, (a, b, a′, b′) = free_associative_algebra(field, ["a", "b", "a'", "b'"])
T = elem_type(algebra)
inverses = Union{T, Tuple{T, T}}[(a, a′), (b, b′)]
relators = Union{T, Tuple{T, Weight}}[ a^9, b^2, (a*b)^4, (a^2*b)^3 ]
subgens = NTuple{1, T}[]

println("vectorenumeration:")
@btime vector_enumeration(algebra, inverses, relators, subgens)

idealgen = [a*a′ - 1, b*b′ - 1, a^9 - 1, b^2 - 1, (a*b)^4 - 1, (a^2*b)^3 - 1]

println("groebnerbasis:")
@btime Generic.groebner_basis(idealgen)

#PSL₃(4)
field = GF(2)
algebra, (a, b, a′, b′) = free_associative_algebra(field, ["a", "b", "a'", "b'"])
T = elem_type(algebra)
inverses = Union{T, Tuple{T, T}}[(a, a′), (b, b′)]
relators = Union{T, Tuple{T, Weight}}[ a^5, b^3, (a*b)^4, (a′*b*a′*b′*a*b)^3, (b*a^2*b*a′^2*b′*a′^2*b′*a^2*b*a*b′*a′) ]
subgens = NTuple{1, T}[]

println("vectorenumeration:")
@btime vector_enumeration(algebra, inverses, relators, subgens)

idealgen = [a*a′ - 1, b*b′ - 1, a^5 - 1, b^3 - 1, (a*b)^4 - 1, (a′*b*a′*b′*a*b)^3 - 1, (b*a^2*b*a′^2*b′*a′^2*b′*a^2*b*a*b′*a′) - 1 ]

println("groebnerbasis:")
@btime Generic.groebner_basis(idealgen)

#M₁₁⁽¹⁾
field = GF(2)
algebra, (a, b, c, a′, b′, c′) = free_associative_algebra(field, ["a", "b", "c", "a'", "b'", "c'"])
T = elem_type(algebra)
inverses = Union{T, Tuple{T, T}}[(a, a′), (b, b′), (c, c′)]
relators = Union{T, Tuple{T, Weight}}[ a^11, b^5, c^4, (a^4*c^2)^3, (b*c^2)^2, (a*b*c)^3, b′*a*b*a′^4, c′*b*c*b′^2 ]
subgens = NTuple{1, T}[]

println("vectorenumeration:")
@btime vector_enumeration(algebra, inverses, relators, subgens)

idealgen = [a*a′ - 1, b*b′ - 1,c*c′ - 1, a^11 - 1, b^5 - 1, c^4 - 1, (a^4*c^2)^3 - 1, (b*c^2)^2 - 1, (a*b*c)^3 - 1, b′*a*b*a′^4 - 1, c′*b*c*b′^2 - 1]

println("groebnerbasis:")
@btime Generic.groebner_basis(idealgen)

#M₁₁⁽²⁾
field = GF(2)
algebra, (a, b, c, a′, b′, c′) = free_associative_algebra(field, ["a", "b", "c", "a'", "b'", "c'"])
T = elem_type(algebra)
inverses = Union{T, Tuple{T, T}}[(a, a′), (b, b′), (c, c′)]
relators = Union{T, Tuple{T, Weight}}[ a^11, b^5, c^4, (a*c)^3, c′*b*c*b′^2, b′*a*b*a′^4 ]
subgens = NTuple{1, T}[]

println("vectorenumeration:")
@btime vector_enumeration(algebra, inverses, relators, subgens)

idealgen = [a*a′ - 1, b*b′ - 1,c*c′ - 1, a^11 - 1, b^5 - 1, c^4 - 1, (a*c)^3 - 1, c′*b*c*b′^2 - 1, b′*a*b*a′^4 - 1 ]

println("groebnerbasis:")
@btime Generic.groebner_basis(idealgen)

#Neu
field = GF(2)
algebra, (a, b, c, a′, b′, c′) = free_associative_algebra(field, ["a", "b", "c", "a'", "b'", "c'"])
T = elem_type(algebra)
inverses = Union{T, Tuple{T, T}}[(a, a′), (b, b′), (c, c′)]
relators = Union{T, Tuple{T, Weight}}[ a^3, b^3, c^3, (a*b)^5, (a′*b)^5, (a*c)^4, (a*c′)^4, a*b′*a*b*c′*a*c*a*c′, (b*c)^3, (b′*c)^4 ]
subgens = NTuple{1, T}[]

println("vectorenumeration:")
@btime vector_enumeration(algebra, inverses, relators, subgens)

idealgen = [a*a′ - 1, b*b′ - 1,c*c′ - 1, a^3 - 1, b^3 - 1, c^3 - 1, (a*b)^5 - 1, (a′*b)^5 - 1, (a*c)^4 - 1, (a*c′)^4 - 1, a*b′*a*b*c′*a*c*a*c′ - 1, (b*c)^3 - 1, (b′*c)^4 - 1 ]

println("groebnerbasis:")
@btime Generic.groebner_basis(idealgen)

#Weyl B₆
field = GF(2)
algebra, (a, b, c, d, e, f) = free_associative_algebra(field, ["a", "b", "c", "d", "e", "f"])
T = elem_type(algebra)
inverses = Union{T, Tuple{T, T}}[a, b, c, d, e, f]
relators = Union{T, Tuple{T, Weight}}[ (a*b)^3, (a*c)^2, (a*d)^2, (a*e)^2, (a*f)^2, (b*c)^3, (b*d)^2, (b*e)^2, (b*f)^2, (c*d)^3, (c*e)^2, (c*f)^2, (d*e)^3, (d*f)^2, (e*f)^4 ]
subgens = NTuple{1, T}[]

println("vectorenumeration:")
@btime vector_enumeration(algebra, inverses, relators, subgens)

idealgen = [a^2 - 1, b^2 - 1, c^2 - 1, d^2 - 1, e^2 - 1, f^2 - 1, (a*b)^3 - 1, (a*c)^2 - 1, (a*d)^2 - 1, (a*e)^2 - 1, (a*f)^2 - 1, (b*c)^3 - 1, (b*d)^2 - 1, (b*e)^2 - 1, (b*f)^2 - 1, (c*d)^3 - 1, (c*e)^2 - 1, (c*f)^2 - 1, (d*e)^3 - 1, (d*f)^2 - 1, (e*f)^4 - 1]

println("groebnerbasis:")
@btime Generic.groebner_basis(idealgen)


