import Base: getindex, setindex!


"""
    getweight(b::BasisElement) 

Return the weight at which the basis element `b` was defined.
"""
getweight(b::BasisElement) = b.weight


"""
    isactiondefined(b::BasisElement, gen::Generator)

Return `true` if the action `b.x` is defined.
"""
isactiondefined(b::BasisElement, x::Generator) = isassigned(b.entries, index(x))


"""
    getindex(b::BasisElement, gen::Generator)

Return the value of the action `b.x` or fail if `b.x` is undefined.
"""
getindex(b::BasisElement, gen::Generator) = b.entries[index(gen)]


"""
    getindex(b::BasisElement, x::Int)

Return the value of the action `b.x` or fail if `b.x` is undefined.
"""
getindex(b::BasisElement, x::Int) = b.entries[x]


"""
    setindex!(b::BasisElement, v::AbstractBasisVector, x::Generator)

Define the action `b.x` as `v`.
"""
function setindex!(b::BasisElement{T}, v::AbstractBasisVector{T}, x::Generator) where {T}
    @assert !isdeleted(b) && !isactiondefined(b, x)  # can't define on deleted basis elements or redefine already defined action

    @debug "defining $b.$x := $v"
    b.entries[index(x)] = v
    b.undefined -= 1
    if b.undefined == 0
        basis(b).completed += 1
    end
end


"""
    define_action!(b::BasisElement, x::Generator, weight::Int, generators::Vector{Generator})::BasisElement

Define the action `b.x` by creating a new basis element `b'` at the weight `weight` and return `b'`. 
"""
function define_action!(b::BasisElement, x::Generator, weight::Int, generators::Vector{Generator})
    @assert !isdeleted(b)

    # create a new basis element b' to define the action 
    img = b.image
    img′ = img === nothing ? nothing : (img[1], img[2] * element(x))
    b′ = create!(basis(b), weight, length(generators), img′)
    @debug "creating basis element $b′ with image $img′ and weight $weight to define $b.$x, resulting in basis size $(nalive(basis(b))) / $(ncreated(basis(b)))"

    # define b.x := b' and if x is invertible with inverse x⁻¹ also b'.x⁻¹ := b to shortcut relations for inverse elements
    # FIXME: for some reason this is not type stable? is this just to be expected with arrays of abstract types?
    b[x] = vector(b′)
    if is_invertible(x)
        b′[inverse(x)] = vector(b)
    end

    return b′
end


# FIXME: How to make this type stable?
struct GeneratorEquation
    # stackable information that v.x = w
    v::AbstractBasisVector
    x::Generator
    w::AbstractBasisVector
end

"""
    replace!(b::BasisElement, replacement::AbstractBasisVector, weight::Union{Int,Missing}, generators::Vector{Generator}, stack)

Delete the basis element `b` and replace it with `replacement` at the weight `weight` and return success.

The replacement must have entries with indices lower than `b`.

Resulting generator equations will not be processed directly but placed on `stack` for later processing.
"""
function replace!(b::BasisElement{T}, replacement::AbstractBasisVector{T}, weight::Union{Int,Missing}, generators::Vector{Generator}, stack) where {T}
    @assert !isdeleted(b)

    @debug "deleting basis element $b and replacing with $replacement, resulting in basis size $(nalive(basis(b))-1) / $(ncreated(basis(b)))"
    delete!(b, replacement)

    # for every x ∈ X with b.x ≠ ⟂ we obtain the generator equation replacement.x = (b.x)
    # these are stacked rather than processed right away to prevent recursion (which may otherwise lead to stack overflow during a massive collapse)
    for x in generators
        if isactiondefined(b, x)
            push!(stack, GeneratorEquation(replacement, x, b[x]))
        end
    end

    # if all entries were filled total number of completed basis elements decrements
    if b.undefined == 0
        basis(b).completed -= 1
    end

    # cleanup
    b.undefined = 0
    empty!(b.entries)

    return true  # success
end


"""
    action!(w::AbstractBasisVector, x::Generator, weight::Union{Int, Missing}, generators::Vector{Generator}, defining::Bool=true)

Compute the action `w.x` with defining if `defining` is set for a generator `x` at the weight `weight` or return `nothing` if this fails.

If `defining` is disabled, no new basis elements will be defined whilst computing `w.x`, possibly leading to failure.
"""
function action!(w::AbstractBasisVector, x::Generator, weight::Union{Int,Missing}, generators::Vector{Generator}, defining::Bool=true)
    result = zero(w)

    for t in w
        b = elem(t)
        λ = coeff(t)
        if isactiondefined(b, x)
            add_scaled!(result, b[x], λ)
        else
            defining || return nothing
            add_scaled!(result, define_action!(b, x, weight, generators), λ)
        end
    end

    return result
end


"""
    action!(w::AbstractBasisVector, g::GroupWord, weight::Union{Int, Missing}, generators::Vector{Generator}, defining::Bool=true)

Compute the action `w.g` with defining if `defining` is set for a group word `g` at the weight `weight` or return `nothing` if this fails.

If `defining` is disabled, no new basis elements will be defined whilst computing `w.g`, possibly leading to failure.
"""
function action!(w::AbstractBasisVector, g::GroupWord, weight::Union{Int,Missing}, generators::Vector{Generator}, defining::Bool=true)
    result = w

    for x in g
        tmp = action!(result, x, weight, generators, defining)
        isnothing(tmp) && return nothing
        result = tmp
    end

    return result
end


"""
    action!(w::AbstractBasisVector, a::AlgebraWord, weight::Union{Int, Missing}, generators::Vector{Generator}, defining::Bool=true)

Compute the action `w.a` with defining if `defining` is set for an algebra word `a` at the weight `weight` or return `nothing` if this fails.

If `defining` is disabled, no new basis elements will be defined whilst computing `w.a`, possibly leading to failure.
"""
function action!(w::AbstractBasisVector, a::AlgebraWord, weight::Union{Int,Missing}, generators::Vector{Generator}, defining::Bool=true)
    result = zero(w)

    for (μ, g) in a
        tmp = action!(w, g, weight, generators, defining)
        isnothing(tmp) && return nothing
        add_scaled!(result, tmp, μ)
    end

    return result
end


"""
    actioncost(v::AbstractBasisVector, x::Generator, weight::Union{Int, Missing})

Calculate the cost of computing the action with defining `v.a` at the weight `weight`.

Returns zero if and only if `v.a` can be computed without defining new basis elements. 
"""
function actioncost(v::AbstractBasisVector, x::Generator, weight::Union{Int,Missing})
    cost = 0
    for t in v
        b = elem(t)
        if !isactiondefined(b, x)
            # there are two options for computing the costs for definitions with weights:
            # (alternatively the number of definitions required could be counted without weights)
            # option 1: compute maximal cost (note: vectors are iterated by index decreasing, hence by weight decreasing)
            cost = getweight(b)
            break
            # option 2: compute total cost
            # cost₁ += getweight(b)
        end
    end

    return cost
end


"""
    process_equation!(v::AbstractBasisVector, weight::Union{Int,Missing}, generators::Vector{Generator}, stack)

Process the information `v = 0` at the weight `weight` whilst transforming `v` in-place and return success.

Resulting generator equations will not be processed directly but placed on `stack` for later processing.
"""
function process_equation!(v::AbstractBasisVector, weight::Union{Int,Missing}, generators::Vector{Generator}, stack)
    @debug "processing equation $v = 0 at weight $weight"

    iszero(v) && return true  # the equation has already been satisfied

    # find the last undeleted basis element b with a non-zero coefficient λ in v to obtain the coincidence b = (-1 / λ) (v - λ⋅b) = b - v / λ =: v, 
    # which can be processed by deleting b and replacing it with v
    t = last(v)
    b = elem(t)
    λ = coeff(t)

    scale!(v, -1 / λ)
    add!(v, b)

    return replace!(b, v, weight, generators, stack)
end


"""
    process_equation!(v₁::AbstractBasisVector, v₂::AbstractBasisVector, weight::Union{Int, Missing}, generators::Vector{Generator})

Process the information `v₁ = v₂` at the weight `weight` and return success.

Resulting generator equations will not be processed directly but placed on `stack` for later processing.
"""
# process_equation!(v₁::AbstractBasisVector, v₂::AbstractBasisVector, weight::Union{Int, Missing}, generators::Vector{Generator}, stack) = process_equation!(v₁ - v₂, weight, generators, stack)


"""
    process_generator_equation!(v₁::AbstractBasisVector, x::Generator, v₂::AbstractBasisVector, weight::Union{Int,Missing}, generators::Vector{Generator}, stack, defining::Bool=true)

Process the information `v₁.x = v₂` at the weight `weight` and return success.

If `defining` is disabled, no new basis elements will be defined whilst computing `v₁.x`, possibly leading to failure.

Further resulting generator equations will be pushed onto `stack` for later processing.
"""
function process_generator_equation!(v₁::AbstractBasisVector{T}, x::Generator, v₂::AbstractBasisVector{T}, weight::Union{Int,Missing}, generators::Vector{Generator}, stack, defining::Bool=true) where {T}
    @debug "processing generator equation $v₁.$x = $v₂ at weight $weight"

    # try to compute v₁.x without defining new basis elements to obtain the equation v₁.x - v₂ = 0
    result = zero(v₁)
    b₀ = nothing
    λ₀ = nothing

    for t in v₁
        b = elem(t)
        λ = coeff(t)

        if isactiondefined(b, x)
            add_scaled!(result, b[x], λ)
        elseif isnothing(b₀)
            # this is the first undefined action encountered, store it but skip defining and adding λ⋅b.x to the result,
            # as from now on v₁.x - λ⋅b.x is computed to obtain the deduction b.x = -(1/λ)(v₁.x - λ⋅b.x - v₂) = -(1/λ)⋅(result - v₂)
            # FIXME: this always gets the highest b, maybe keep on updating to get the lowest b for more valuable information?
            b₀ = b
            λ₀ = λ
        else
            defining || return false
            add_scaled!(result, define_action!(b, x, weight, generators), λ)
        end
    end

    if isnothing(b₀)
        # process the equation 0 = v₁.x - v₂
        subtract!(result, v₂)
        return process_equation!(result, weight, generators, stack)
    else
        # process the deduction b₀.x = -(1/λ₀)(v₁.x - λ₀⋅b₀.x - v₂)
        subtract!(result, v₂)
        scale!(result, -1 / λ₀)

        # define b₀.x = result and if x is invertible with inverse x⁻¹ also stack/process the information that result.x⁻¹ = b₀ to make the relations xx⁻¹ = 1 and x⁻¹x = 1 redundant
        b₀[x] = result
        if is_invertible(x)
            push!(stack, GeneratorEquation(result, inverse(x), vector(b₀)))
        end

        return true  # success
    end
end


"""
    process_subgen!(subgen::NTuple{n,Union{AlgebraWord,GroupWord,Generator}}, start_basis::Vector{BasisElement{T}}, weight::Int, generators::Vector{Generator}) where {T, n}

Process the information `b₁.s₁ + ⋯ + bₙ.sₙ = 0` for the submodule generator `subgen = (s₁, …, sₙ)` and the initial basis elements `b₁, …, bₙ` in `start_basis` at the weight `weight`.
"""
function process_subgen!(subgen::NTuple{n,Union{AlgebraWord,GroupWord,Generator}}, start_basis::Vector{BasisElement{T}}, weight::Int, generators::Vector{Generator}) where {T,n}
    @assert length(start_basis) == n  # dimensions must match

    @debug "processing submodule generator $subgen at weight $weight"

    # compute r = b₁.s₁ + ⋯ + bₙ.sₙ
    r = zero(start_basis[1])
    for (b, s) in zip(start_basis, subgen)
        add!(r, action!(vector(b), s, weight, generators, true))
    end

    stack = Vector{GeneratorEquation}()

    # process the equation r = 0
    process_equation!(r, weight, generators, stack)

    # process stacked generator equations
    while !isempty(stack)
        ge = pop!(stack)
        process_generator_equation!(ge.v, ge.x, ge.w, weight, generators, stack, true)
    end
end


### process relations
"""
    process_relation!(b::BasisElement, relation::AlgebraRelation, weight::Union{Int,Missing}, generators::Vector{Generator}, stack, defining::Bool=true)

Process the information `b.a = b` for an algebra word `a` given by the algebra relation `relation` at the weight `weight` and return success.

If `defining` is disabled, no new basis elements will be defined whilst computing `b.a`, possibly leading to failure.

Resulting generator equations may be pushed onto `stack` for later processing.
"""
function process_relation!(b::BasisElement, relation::AlgebraRelation, weight::Union{Int,Missing}, generators::Vector{Generator}, stack, defining::Bool=true)
    @assert !isdeleted(b)  # relations should not be applied to already deleted basis elements

    # most primitive way to apply relators, when more information is available on the relator, use a better method
    a = relation.word
    @debug "processing fixing relation $b.$a = $b at weight $weight"

    # compute b.a, allow defining based on mode
    result = action!(vector(b), a, weight, generators, defining)
    isnothing(result) && return false

    # process the equation b.a = b ⟺ b.a - b = 0
    subtract!(result, b)
    return process_equation!(result, weight, generators, stack)
end


"""
    process_relation!(b::BasisElement, relation::LinearRelation, weight::Union{Int,Missing}, generators::Vector{Generator}, stack, defining::Bool=true)

Process the information `b.a = 0` for a linear algebra word `a` given by the linear relation `relation` at the weight `weight` and return success.

If `defining` is disabled, no new basis elements will be defined whilst computing `b.a`, possibly leading to failure.

Resulting generator equations may be pushed onto `stack` for later processing.
"""
function process_relation!(b::BasisElement, relation::LinearRelation, weight::Union{Int,Missing}, generators::Vector{Generator}, stack, defining::Bool=true)
    @assert !isdeleted(b)  # relations should not be applied to already deleted basis elements

    a = relation.word
    @debug "processing linear relation $b.$a = 0 at weight $weight"

    # try to compute b.a = b.(λ₁x₁ + ⋯ + λₗxₗ + λₗ₊₁) to generate the equation 0 = b.a, if this fails at some xᵢ because b.xᵢ is undefined
    # instead produce the deduction b.xᵢ = -1/λᵢ b.(λ₁x₁ + ⋯ + λᵢ₋₁xᵢ₋₁ + λᵢ₊₁xᵢ₊₁ + ⋯ + λₗxₗ + λₗ₊₁)
    result = zero(b)
    x₀ = nothing
    λ₀ = nothing

    for (λ, g) in a
        if isone(g)
            add_scaled!(result, b, λ)
        else
            x = only(g)  # must be linear

            if isactiondefined(b, x)
                add_scaled!(result, b[x], λ)
            elseif isnothing(x₀)
                # this is the first undefined action encountered, store it but skip defining and adding λ⋅b.x to the result,
                # as from now on b.a - b.(λ⋅x) is computed to obtain the deduction b.x = -(1/λ)(b.a - b.(λ⋅x)) = -(1/λ)⋅result
                x₀ = x
                λ₀ = λ
            else
                defining || return false
                add_scaled!(result, define_action!(b, x, weight, generators), λ)
            end
        end
    end

    if isnothing(x₀)
        return process_equation!(result, weight, generators, stack)
    else
        scale!(result, -1 / λ₀)
        b[x₀] = result

        if is_invertible(x₀)
            push!(stack, GeneratorEquation(result, inverse(x₀), b))
        end

        return true  # success
    end
end


"""
    process_relation!(b::BasisElement, relation::BinomialRelation, weight::Union{Int,Missing}, generators::Vector{Generator}, stack, defining::Bool=true)

Process the information `b.(λ₁g₁) = b.(λ₂g₂)` for two group words `g₁, g₂` with coefficents `λ₁, λ₂` given by the binomial relation `relation` at the weight `weight` and return success.

If `defining` is disabled, no new basis elements will be defined whilst computing `b.(λ₁g₁)` and `b.(λ₂g₂)`, possibly leading to failure. 

Resulting generator equations may be pushed onto `stack` for later processing.
"""
function process_relation!(b::BasisElement, relation::BinomialRelation, weight::Union{Int,Missing}, generators::Vector{Generator}, stack, defining::Bool=true)
    @assert !isdeleted(b)  # relations should not be applied to already deleted basis elements

    λ₁ = relation.coeff₁
    λ₂ = relation.coeff₂
    g₁ = relation.word₁
    g₂ = relation.word₂
    @debug "processing binomial relation ($λ₁⋅$b).$g₁ = ($λ₂⋅$b).$g₂ at weight $weight"

    l₁ = length(g₁)
    l₂ = length(g₂)

    # TODO: replace scalar * basis element with init of SparseVector
    if l₁ == 0 && l₂ == 0
        # process equation λ₁⋅b = λ₂⋅b ⟺ result := λ₁⋅b + (-λ₂)⋅b = 0 
        result = zero(b)
        add_scaled!(result, b, λ₁)
        add_scaled!(result, b, -λ₂)
        return process_equation!(result, weight, generators, stack)
    end

    # result₁ := λ₁ * b
    result₁ = zero(b)
    add_scaled!(result₁, b, λ₁)

    # result₂ := λ₁ * b
    result₂ = zero(b)
    add_scaled!(result₂, b, λ₂)

    l = l₁ + l₂ + 1
    i₁ = 1
    i₂ = 1
    # note: upon reaching the opposite word the inverse of a non-invertible generator and thus one xᵢ can be nothing, 
    # leading to an action cost costᵢ of missing (which is treated like ∞)
    # then the computation will force its way through the remainders opposite side, where no inverses will be taken
    x₁ = i₁ ≤ l₁ ? g₁[i₁] : inverse(g₂[l-i₁])  # next generator x₁ to apply to result₁
    x₂ = i₂ ≤ l₂ ? g₂[i₂] : inverse(g₁[l-i₂])  # next generator x₂ to apply to result₂
    cost₁ = isnothing(x₁) ? missing : actioncost(result₁, x₁, weight)  # definition costs of applying x₁ to result₁
    cost₂ = isnothing(x₂) ? missing : actioncost(result₂, x₂, weight)  # definition costs of applying x₂ to result₂

    # push generators from both sides until both ends meet on the same (last) generator
    while i₁ + i₂ < l
        if isless(cost₁, cost₂)
            # less expensive to push generator on result₁
            defining || cost₁ === 0 || return false  # if no new basis elements may be defined, the cost must be 0 to proceed
            ismissing(cost₁) && return false  # cost must be finite to proceed (otherwise x may be nothing) 
            result₁ = action!(result₁, x₁, weight, generators, defining)  # set result₁ := result₁.x₁
            i₁ += 1
            x₁ = i₁ ≤ l₁ ? g₁[i₁] : inverse(g₂[l-i₁])  # next generator x₁ to apply to result₁
        else
            # less or equal as expensive to push generator on result₂
            defining || cost₂ === 0 || return false  # if no new basis elements may be defined, the cost must be 0 to proceed
            ismissing(cost₂) && return false  # cost must be finite to proceed (otherwise x may be nothing)
            result₂ = action!(result₂, x₂, weight, generators, defining)  # set result₂ := result₂.x₂
            i₂ += 1
            x₂ = i₂ ≤ l₂ ? g₂[i₂] : inverse(g₁[l-i₂])  # next generator x₂ to apply to result₂
        end

        # recompute costs
        cost₁ = isnothing(x₁) ? missing : actioncost(result₁, x₁, weight)  # definition costs of applying x₁ to result₁
        cost₂ = isnothing(x₂) ? missing : actioncost(result₂, x₂, weight)  # definition costs of applying x₂ to result₂
    end

    # try to push last remaining generator x₁ or x₂
    if cost₁ === 0
        # can push generator x₁ at no cost on result₁, i.e. without defining, resulting in equation (result₁.x₁) = result₂ ⟺ (result₁.x₁) - result₂ = 0
        r = action!(result₁, x₁, weight, generators, false)
        subtract!(r, result₂)
        return process_equation!(r, weight, generators, stack)
    elseif cost₂ === 0
        # can push generator x₂ at no cost on result₂, i.e. without defining, resulting in equation (result₂.x₂) = result₁ ⟺ (result₂.x₂) - result₁ = 0
        r = action!(result₂, x₂, weight, generators, false)
        subtract!(r, result₁)
        return process_equation!(r, weight, generators, stack)
    else
        # last generator can't be pushed at no cost, i.e. without defining ...
        defining || return false

        if isless(cost₁, cost₂)
            # less expensive to push generator x₁ on result₁, resulting in generator equation result₁.x₁ = result₂
            return process_generator_equation!(result₁, x₁, result₂, weight, generators, stack, true)
        else
            # less or equal as expensive to push generator x₂ on result₂, resulting in generator equation result₂.x₂ = result₁
            return process_generator_equation!(result₂, x₂, result₁, weight, generators, stack, true)
        end
    end
end


"""
    process_relation!(b::BasisElement, relation::GroupTypeRelation, weight::Union{Int,Missing}, generators::Vector{Generator}, stack, defining::Bool=true)

Process the information `b.g = b` for a group word `g` given by the group type relation `relation` at the weight `weight` and return success.

If `defining` is disabled, no new basis elements will be defined whilst computing `b.g`, possibly leading to failure.

Resulting generator equations may be pushed onto `stack` for later processing.
"""
function process_relation!(b::BasisElement, relation::GroupTypeRelation, weight::Union{Int,Missing}, generators::Vector{Generator}, stack, defining::Bool=true)
    @assert !isdeleted(b)  # relations should not be applied to already deleted basis elements

    g = relation.word
    @debug "processing fixing group-type relation $b.$g = $b at weight $weight"

    l = length(g)
    l == 0 && return true  # relation is trivially satisfied

    i₁ = 1
    i₂ = l
    result₁ = vector(b)
    result₂ = vector(b)
    x₁ = g[i₁]  # next generator x₁ to apply to result₁
    x₂ = inverse(g[i₂])  # next generator x₂ to apply to result₂
    cost₁ = actioncost(result₁, x₁, weight)  # definition costs of applying x₁ to result₁
    cost₂ = isnothing(x₂) ? missing : actioncost(result₂, x₂, weight)  # definition costs of applying x₂ to result₂

    # push generators from both sides until both ends meet on the same (last) generator
    while i₁ < i₂
        if isless(cost₂, cost₁)
            # less expensive to push generator x₂ on result₂
            defining || cost₂ === 0 || return false  # if no new basis elements may be defined, the cost must be 0 to proceed

            result₂ = action!(result₂, x₂, weight, generators, defining)
            i₂ -= 1
            x₂ = inverse(g[i₂])  # next generator x₂ to apply to result₂
        else
            # less or equal as expensive to push generator x₁ on result₁
            defining || cost₁ === 0 || return false  # if no new basis elements may be defined, the cost must be 0 to proceed

            result₁ = action!(result₁, x₁, weight, generators, defining)
            i₁ += 1
            x₁ = g[i₁]  # next generator x₁ to apply to result₁
        end

        # recompute costs
        cost₁ = actioncost(result₁, x₁, weight)  # definition costs of applying x₁ to result₁
        cost₂ = isnothing(x₂) ? missing : actioncost(result₂, x₂, weight)  # definition costs of applying x₂ to result₂
    end

    # try to push last remaining generator x₁ or x₂
    if cost₁ === 0
        # can push generator x₁ at no cost on result₁, i.e. without defining, resulting in equation (result₁.x₁) = result₂ ⟺ (result₁.x₁) - result₂ = 0
        r = action!(result₁, x₁, weight, generators, false)
        subtract!(r, result₂)
        return process_equation!(r, weight, generators, stack)
    elseif cost₂ === 0
        # can push generator x₂ at no cost on result₂, i.e. without defining, resulting in equation (result₂.x₂) = result₁ ⟺ (result₂.x₂) - result₁ = 0
        r = action!(result₂, x₂, weight, generators, false)
        subtract!(r, result₁)
        return process_equation!(r, weight, generators, stack)
    else
        # last generator can't be pushed at no cost, i.e. without defining ...
        defining || return false

        if isless(cost₂, cost₁)
            # less expensive to push generator x₂ on result₂, resulting in generator equation result₂.x₂ = result₁
            return process_generator_equation!(result₂, x₂, result₁, weight, generators, stack, true)
        else
            # less or equal as expensive to push generator x₁ on result₁, resulting in generator equation result₁.x₁ = result₂
            return process_generator_equation!(result₁, x₁, result₂, weight, generators, stack, true)
        end
    end
end


"""
    process_relation!(b::BasisElement, relation::DefineRelation, weight::Union{Int,Missing}, generators::Vector{Generator}, stack, defining::Bool=true)

Verify the action `b.x` is defined for the generator `x` given by the define relation `relation` at the weight `weight` and return success.

If `defining` is disabled, no new basis element will be defined to define `b.x`, possibly leading to failure. 

As this relation type does not result in generator equations, `stack` remains unused.
"""
function process_relation!(b::BasisElement, relation::DefineRelation, weight::Union{Int,Missing}, generators::Vector{Generator}, stack, defining::Bool=true)
    @assert !isdeleted(b)  # relations should not be applied to already deleted basis elements

    x = relation.x
    @debug "verifying $b.$x is defined at weight $weight"

    isactiondefined(b, x) && return true  # relation already satisfied
    defining || return false
    define_action!(b, x, weight, generators)
    return true  # success
end


# TODO: cleanup and document, maybe also include lookahead here, otherwise add functions in another location
function process_relation!(basis::Basis{BasisElement{T}}, relation::Relation, weight::Union{Int,Missing}, generators::Vector{Generator}) where {T}
    @assert isless(0, weight)  # weight must be positive

    # for the relation r process all basis elements b such that getweight(b) + getweight(r) ≤ weight, i.e. getweight(b) ≤ weight - getweight(r) =: blimit 
    # note: missing is treated like ∞
    blimit = weight - getweight(relation)

    # start at the first or continue at the next unprocessed basis element 
    b = isnothing(relation.processed_last) ? basis.first_undeleted : next(relation.processed_last)

    stack = Vector{GeneratorEquation}()

    while !isnothing(b) && isless(getweight(b), blimit)
        process_relation!(b, relation, weight, generators, stack, true) || error("something went wrong")

        # process stacked generator equations; here defining is always okay, even in lookahead
        while !isempty(stack)
            ge::GeneratorEquation = pop!(stack)
            process_generator_equation!(ge.v, ge.x, ge.w, weight, generators, stack, true)
        end

        relation.processed_last = b
        b = next(b)
    end
end
