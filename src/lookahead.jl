import Base: <=, popfirst!, push!, /

function action_lookahead!(w::AbstractBasisVector, x::Generator)
    result = zero(w)

    for (b, λ) in w
        if isdefined(b, x)
            add_scaled!(result, b[x], λ)
        else
            return nothing
        end
    end
    return result
end


function action_lookahead!(w::AbstractBasisVector, g::GroupWord)
    result = w
    for x in g
        result = action_lookahead!(result, x)
        if isnothing(result)
            return nothing
        end
    end
    return result
end

function action_lookahead!(w::AbstractBasisVector, a::AlgebraWord)
    result = zero(w)

    for (μ, g) in a
        tmp = action_lookahead!(w, g)
        if isnothing(tmp)
            return nothing
        end
        add_scaled!(result, tmp, μ)
    end
    return result
end

function process_relation_lookahead!(b::BasisElement, relator::AlgebraRelation, weight::Union{Int,Missing}; generators::Vector{Generator})
    @assert !isdeleted(b)  # relations should not be applied to already deleted basis elements
    a = relator.word
    @debug "trying to process fixing relation $b.$a = $b in lookahead-mode"

    res = action_lookahead!(b, a)
    if !isnothing(res)
        @debug "able to process in lookahead-mode"
        process_equation!(res, b, weight; generators)
        return true
    end
    @debug "unable to process in lookahead-mode"
    return false
end

function process_relation_lookahead!(b::BasisElement, relator::GroupTypeRelation, weight::Union{Int,Missing}; generators::Vector{Generator})
    @assert !isdeleted(b)
    a = relator.word
    @debug "trying to process fixing relation $b.$a = $b in lookahead-mode"

    i₁ = 0
    i₂ = length(a.word)

    v₁ = v₂ = b

    while i₁ <= i₂
        tmp = action_lookahead!(v₁, a[i₁])
        if !isnothing(tmp)
            i₁ += 1
            v₁ = tmp
        elseif is_invertible(a[i₂])
            tmp = action_lookahead!(v₂, inverse(a[i₂]))
            if isnothing(tmp)
                @debug "unable to process in lookahead-mode"
                return false
            else
                v₂ = tmp
                i₂ -= 1
            end
        else
            @debug "unable to process in lookahead-mode"
            return false
        end
    end

    process_equation!(v₁, v₂, weight; generators)
    return true
end

function process_relation_lookahead!(b::BasisElement, relator::BinomialRelation, weight::Union{Int,Missing}; generators::Vector{Generator})
    @assert !isdeleted(b)  # relations should not be applied to already deleted basis elements

    λ₁ = relator.coeff₁
    λ₂ = relator.coeff₂
    g₁ = relator.word₁
    g₂ = relator.word₂
    @debug "trying to binomial relation ($λ₁⋅$b).$g₁ = ($λ₂⋅$b).$g₂ at weight $weight"

    result₁ = λ₁ * b
    result₂ = λ₂ * b

    l₁ = length(g₁)
    l₂ = length(g₂)
    @assert l₁ > 0 || l₂ > 0  # not both words can be one
    l = l₁ + l₂ + 1
    i₁ = 1
    i₂ = 1

    x₁ = i₁ <= l₁ ? g₁[i₁] : inverse(g₂[l-i₁])
    x₂ = i₂ <= l₂ ? g₂[i₂] : inverse(g₁[l-i₂])
    tmp₁ = action_lookahead!(result₁, x₁)
    tmp₂ = action_lookahead!(result₂, x₂)
    while isless(i₁ + i₂, l)
        if !isnothing(tmp₁)
            i₁ += 1
            x₁ = i₁ <= l₁ ? g₁[i₁] : inverse(g₂[l-i₁])
            result₁ = tmp₁
            tmp₁ = action_lookahead!(result₁, x₁)
        elseif !isnothing(tmp₂)
            i₂ += 1
            x₂ = i₂ <= l₂ ? g₂[i₂] : inverse(g₁[l-i₂])
            result₂ = tmp₂
            tmp₂ = action_lookahead!(result₂, x₂)
        else
            @debug "unable to process in lookahead-mode"
            return false
        end
    end
    if !isnothing(tmp₁)
        process_equation!(tmp₁, result₂, weight; generators)
        @debug "able to process in lookahead-mode"
        return true
    elseif !isnothing(tmp₂)
        process_equation!(result₁, tmp₂, weight; generators)
        @debug "able to process in lookahead-mode"
        return true
    end
    @debug "unable to process in lookahead-mode"
    return false
end

function process_relation_lookahead!(b::BasisElement, relator::LinearRelation, weight::Union{Int,Missing}; generators::Vector{Generator})
    @assert !isdeleted(b)

    a = relator.word
    @debug "trying to process linear relation $b.$a = 0 at weight $weight in lookahead-mode"

    result = zero(b)
    xᵤ = nothing
    λᵤ = nothing
    no_define = true

    for (μ, g) in a
        if isone(g)
            add_scaled!(result, b, μ)
        else
            x = only(g)
            if !isdefined(b, x)
                if no_define
                    xᵤ = x
                    λᵤ = μ
                    no_define = false
                else
                    @debug "unable to process in lookahead-mode"
                    return false
                end
            else
                add_scaled!(result, b[x], μ)
            end
        end
    end

    @debug "able to process in lookahead-mode"
    if no_define
        process_equation!(result, zero(b), weight; generators)
    else
        define_action!(b, xᵤ, result / -λᵤ, weight; generators)
    end
    return true
end

function process_relation_lookahead!(b::BasisElement, relator::DefineRelation, weight::Union{Int,Missing}; generators::Vector{Generator})
    @debug "not processing DefineRelations in lookahead-mode"
    return false
end




function process_relation_lookahead!(basis::Basis, relation::Relation, weight::Union{Int,Missing}; generators::Vector{Generator})
    @assert isless(0, weight)

    blimit = weight - getweight(relation, lookahead=true)
    #    non_processed = false

    b = isnothing(relation.processed_last) ? basis.first_undeleted : next(relation.processed_last)

    while !isnothing(b) && isless(getweight(b), blimit)
        processed = process_relation_lookahead!(b, relation, weight; generators)
        #Methode: saving missed basis
        #=        
                if !processed
                    push!(relation.missed, b)   #TODO add missed to relation
                    non_processed = true        #? better methode ?
                end
                relation.processed_last = b
        =#
        b = next(b)
    end
    #    return non_processed
end

function process_missed(basis::Basis, relation::Relation, weight::Union{Int,Missing}; generators::Vector{Generator})

    Q = relator.missed
    while !isempty(Q)
        b = popfirst!(Q)
        if !isdeleted(b)
            process_relation(b, relation, weight; generators)
        end
    end

end

#TODO find fitting criterion for lookahead
function activate_lookahead(basis::Basis)
    x = 1 - basis.deleted / basis.created
    if 0 <= x
        return true
    end
    return false
end



