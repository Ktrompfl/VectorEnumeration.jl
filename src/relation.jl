
struct Weight
    normal::Int  # weight in normal mode
    lookahead::Int  # weight in lookahead mode
    function Weight(normal::Int, lookahead::Int)
        if !(normal > 0 && lookahead > 0)
            error("weights must be positive")
        end
        return new(normal, lookahead)
    end
end

Weight(value::Int) = Weight(value, value)


abstract type Relation end

mutable struct AlgebraRelation <: Relation
    # relation word = 1 for a word of the algebra
    const word::AlgebraWord
    const weight::Weight
    processed_last::Union{Nothing,BasisElement}
    AlgebraRelation(word::AlgebraWord, weight::Weight) = new(word, weight, nothing)
end

mutable struct GroupTypeRelation <: Relation
    # relation word = 1 for a group word
    const word::GroupWord
    const weight::Weight
    processed_last::Union{Nothing,BasisElement}
    GroupTypeRelation(word::GroupWord, weight::Weight) = new(word, weight, nothing)
end

mutable struct LinearRelation <: Relation
    # relation word = 0 for a linear (of degree 1) algebra word
    const word::AlgebraWord
    const weight::Weight
    processed_last::Union{Nothing,BasisElement}
    LinearRelation(word::AlgebraWord, weight::Weight) = new(word, weight, nothing)
end

mutable struct BinomialRelation <: Relation
    # relation coeff₁⋅word₁ = coeff₂⋅word₂
    const coeff₁::FieldElement
    const word₁::GroupWord
    const coeff₂::FieldElement
    const word₂::GroupWord
    const weight::Weight
    processed_last::Union{Nothing,BasisElement}
    BinomialRelation(coeff₁::FieldElement, word₁::GroupWord, coeff₂::FieldElement, word₂::GroupWord, weight::Weight) = new(coeff₁, word₁, coeff₂, word₂, weight, nothing)
end

BinomialRelation(word₁::GroupWord, word₂::GroupWord, weight::Weight) = BinomialRelation(one(base_ring(word₁)), word₁, one(base_ring(word₂)), word₂, weight)  # word₁ = word₂
BinomialRelation(word₁::GroupWord, weight::Weight) = BinomialRelation(word₁, GroupWord([], one(parent(word₁))), weight)  # word₁ = 1

mutable struct DefineRelation <: Relation
    # virtual relation x = x for x ∈ X, required to ensure b.x is defined for all b ∈ B and x ∈ X
    const x::Generator
    const weight::Weight
    processed_last::Union{Nothing,BasisElement}
    DefineRelation(x::Generator, weight::Weight) = new(x, weight, nothing)
end


getweight(r::Relation; lookahead::Bool=false) = lookahead ? r.weight.lookahead : r.weight.normal  # TODO: how to correctly do this with abstract types


## parse the relator to obtain the most efficient relation implementation
function Relation(relator::FreeAssAlgElem, gens::Vector{Generator}; weight::Union{Weight,Nothing}=nothing, fixing::Bool=true)
    A = parent(relator)
    R = base_ring(relator)

    # check if relator describes fixing relation relator = 1 or annihilating relation relator = 0 and
    # convert annihilating relation to fixing relation relator = 0 ⟺ r := relator + 1 = 1
    r = fixing ? relator : relator + one(A)

    iszero(r) && error("relation 0 = 1 invalid for algebra over $R")

    if total_degree(r) <= 1
        # linear relation, default weight 3
        a = r - one(A)  # convert r = 1 to annihilating relation a := r - 1 = 0, which is still linear
        word = AlgebraWord(a, gens)
        weight = isnothing(weight) ? Weight(3) : weight

        return LinearRelation(word, weight)
    end

    if length(r) == 1
        word = GroupWord(monomial(r, 1), gens)

        #  (scaled) group type relation λ⋅g := r = 1 default weight length(g) ÷ 2
        weight = isnothing(weight) ? Weight(length(word) ÷ 2) : weight

        # TODO: check for λ = 1 and use GroupTypeRelator or keep this more general?

        return BinomialRelation(coeff(r, 1), word, one(R), GroupWord(one(A), gens), weight)
    elseif length(r) == 2
        word₁ = GroupWord(monomial(r, 1), gens)
        word₂ = GroupWord(monomial(r, 2), gens)

        if isone(word₁)
            # binomial relation / (scaled) group type relation λ₁⋅1 + λ₂⋅g₂ := r = 1 ⟺ λ₂⋅g₂ = (1 - λ₁)⋅1; default weight length(g₂) ÷ 2
            weight = isnothing(weight) ? Weight(length(word₂) ÷ 2) : weight

            return BinomialRelation(coeff(r, 2), word₂, 1 - coeff(r, 1), word₁, weight)
        elseif isone(word₂)
            # binomial relation / (scaled) group type relation λ₁⋅g₁ + λ₂⋅1 := r = 1 ⟺ λ₁⋅g₁ = (1 - λ₂)⋅1; default weight length(g₁) ÷ 2
            weight = isnothing(weight) ? Weight(length(word₁) ÷ 2) : weight

            return BinomialRelation(coeff(r, 1), word₁, 1 - coeff(r, 2), word₂, weight)
        end
    elseif length(r) == 3
        a = r - one(A)  # convert r = 1 to annihilating relation a := r - 1 = 0 and check if a is binomial
        if length(a) == 2
            word₁ = GroupWord(monomial(r, 1), gens)
            word₂ = GroupWord(monomial(r, 2), gens)

            # binomial relation λ₁⋅g₁ + λ₂⋅g₂ := a = 0 ⟺ λ₁⋅g₁ = -λ₂⋅g₂; default weight (length(g₁) + length(g₂)) ÷ 2
            weight = isnothing(weight) ? Weight((length(word₁) + length(word₂)) ÷ 2) : weight

            return BinomialRelation(coeff(r, 1), word₁, -coeff(r, 2), word₂, weight)
        end
    end

    # no special case can be applied, hence default to algebra relation r = 1; default weight 3
    word = AlgebraWord(r, gens)
    weight = isnothing(weight) ? Weight(3) : weight

    return AlgebraRelation(word, weight)
end
