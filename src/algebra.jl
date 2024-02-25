import AbstractAlgebra: coeff, monomial, base_ring, parent, is_invertible
import Base: +, *, -, getindex, length, iszero, isone, eltype, iterate, firstindex, lastindex, show

# wrap different kinds of free associative algebra elements, so they can be matched on and can have extra data attached
mutable struct Generator
    const index::Int
    const element::FreeAssAlgElem
    inverse::Union{Generator,Nothing}
    Generator(index::Int, element::FreeAssAlgElem) = new(index, element, nothing)
end

generators(A::FreeAssAlgebra) = [Generator(i, gen(A, i)) for i in 1:nvars(A)]
isgenerator(a::FreeAssAlgElem) = length(a) == 1 && isone(coeff(a, 1)) && length(exponent_word(a, 1)) == 1

is_invertible(x::Generator) = !isnothing(x.inverse)
inverse(x::Generator) = x.inverse

function Generator(a::FreeAssAlgElem)
    @assert isgenerator(a)
    return Generator(a.exps[1][1], a)
end

struct GroupWord
    word::Vector{Generator}
    element::FreeAssAlgElem
end

isgroupword(a::FreeAssAlgElem) = length(a) == 1 && isone(coeff(a, 1))

function GroupWord(a::FreeAssAlgElem, gens::Vector{Generator})
    @assert isgroupword(a)
    return GroupWord([gens[i] for i in a.exps[1]], a)
end

struct AlgebraWord{R<:FieldElement}
    coeffcients::Vector{R}
    monomials::Vector{GroupWord}
    element::FreeAssAlgElem
end

function AlgebraWord(a::FreeAssAlgElem, gens::Vector{Generator})
    return AlgebraWord(a.coeffs, [GroupWord(m, gens) for m in monomials(a)], a)
end


# only display the element
show(io::IO, x::Generator) = show(io, x.element)
show(io::IO, g::GroupWord) = show(io, g.element)
show(io::IO, a::AlgebraWord) = show(io, a.element)

element(x::Generator) = x.element
element(g::GroupWord) = g.element
element(a::AlgebraWord) = a.element
parent(x::Generator) = parent(x.element)
parent(g::GroupWord) = parent(g.element)
parent(a::AlgebraWord) = parent(a.element)
base_ring(x::Generator) = base_ring(x.element)
base_ring(g::GroupWord) = base_ring(g.element)
base_ring(a::AlgebraWord) = base_ring(a.element)

Base.getindex(g::GroupWord, i::Int) = g.word[i]
Base.firstindex(g::GroupWord) = firstindex(g.word)
Base.lastindex(g::GroupWord) = lastindex(g.word)
Base.length(g::GroupWord) = length(g.word)
Base.eltype(::Type{GroupWord}) = Generator
Base.iterate(g::GroupWord, state::Int=1) = (1 <= state <= length(g)) ? (g[state], state + 1) : nothing
Base.isone(g::GroupWord) = length(g) == 0

Base.getindex(a::AlgebraWord, i::Int) = (a.coeffcients[i], a.monomials[i])
Base.firstindex(a::AlgebraWord) = (firstindex(a.coeffcients), firstindex(a.monomials))
Base.lastindex(a::AlgebraWord) = (lastindex(a.coeffcients), lastindex(a.monomials))
Base.length(a::AlgebraWord) = length(a.coeffcients)
Base.eltype(::Type{AlgebraWord}) = Tuple{FieldElement,GroupWord}
Base.iterate(a::AlgebraWord, state::Int=1) = (1 <= state <= length(a)) ? ((a.coeffcients[state], a.monomials[state]), state + 1) : nothing
Base.isone(a::AlgebraWord) = length(a) == 1 && isone(a.coeffcients[1]) && isone(a.monomials[i])
Base.iszero(a::AlgebraWord) = length(a) == 0

word(g::GroupWord) = g.word
coeff(a::AlgebraWord, i::Int) = a.coeffcients[i]
monomial(a::AlgebraWord, i::Int) = a.monomials[i]
index(g::Generator) = g.index


*(x₁::Generator, x₂::Generator) = GroupWord([x₁, x₂], x₁.element * x₂.element)
*(x::Generator, g::GroupWord) = GroupWord([x, g.word...], x.element * g.element)
*(g::GroupWord, x::Generator) = GroupWord([g.word..., x], g.element * x.element)
*(g₁::GroupWord, g₂::GroupWord) = GroupWord([g₁..., g₂...], g₁.element * g₂.element)

*(λ::R, x::Generator) where {R<:FieldElement} = AlgebraWord{R}([λ], [GroupWord(x, x.element)], λ * x.element)
*(x::Generator, λ::R) where {R<:FieldElement} = AlgebraWord{R}([λ], [GroupWord(x, x.element)], λ * x.element)
*(λ::R, g::GroupWord) where {R<:FieldElement} = AlgebraWord{R}([λ], [g], λ * g.element)
*(g::GroupWord, λ::R) where {R<:FieldElement} = AlgebraWord{R}([λ], [g], g.element * λ)
