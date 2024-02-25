module VectorEnumeration

using AbstractAlgebra  # AbstractAlgebra.jl covers everything needed for the implementation
using SparseArrays  # SparseArrays are used for the output

export vector_enumeration, base_qa, base_qm, dimension_qa, dimension_qm, matrices_qa, matrices_qm, Weight

include("algebra.jl")
include("vector.jl")
include("relation.jl")
include("process.jl")
#include("lookahead.jl") #! outcommented part for lookahead


"""
    vector_enumeration(algebra::Any, inverses::AbstractArray, relators::AbstractArray, rank::Int=1, submodule_generators::AbstractArray=[]; abelian::Bool = false, max_weight::Int = 100)

Compute a matrix representation for a finitely presented algebra over a field given by generators and relations.

# Arguments
- `algebra::FreeAssAlgebra`: the free associative algebra (over a field) in which computations take place.
- `inverses::Vector{Union{A, Tuple{A, A}}}`: the inverses for generators of the algebra. Involutionary generators should be given as single values, other invertible elements as tuple together with their inverse.
    This will automatically adjoin the relations `x^2 = 1` resp. `xx′ = 1` and `x′x = 1` at weight 3 and weight 6 in lookahead.
- `relators::Vector{Union{A, Tuple{A, Weight}}}`: the relations of the algebra as fixing relators.
    Annihilating relators `r` (such that `r = 0`) need to be converted to the respective fixing relators `r′ := r + 1` (such that `r′ = 1`).
    Optionally relators can be grouped together in a pair with a weight, otherwise by default weights are assumed.
    For group-type relators the default weight is half the length of the relator. For an algebra relation, the default weight is 3.
    To ensure the result is complete, for any generator `x` the relation `x = x` is automatically adjoined at weight 3 and weight 6 in lookahead mode.
- `rank::Int`: the amount of generators of the finitely presented module.
- `submodule_generators::Vector{NTuple{S, FreeAssAlgElem}}`: a list of submodule generators.
- `abelian::Bool = false`: Indicate that all generators may be assumed to commute with each other, e.g. the algebra is abelian, hence a (multivariate) polynomial ring, to automatically adjoin the relation `xy = yx` at weight 2 for every pair of generators.
- `max_weight::Int = 100`: The maximum weight used during computation. Upon exceeding this weight without the calculation being complete, the program fails. 

# Output
- `variables`: the generating variables of `algebra`. \n
- `matrices` : the images of the `variables` under the computed representation in order of `variables` as sparse matrices. \n
- `dimension`: the dimension of the computed representation. \n
- `ntotal`   : the total number of created basis elements. \n
- `images`   : \n
- `preimages`: \n

# Example
Considering the algebra `A = < a, b, c | a², b², c², (ab)³, (bc)³, (ac)²>ₖ ≅ kS₃` with `k = F₇`
and the module `M = < e | e.abc = e >` over `A`

Input:

`algebra, (a, b, c) = free_associative_algebra(GF(7), ["a", "b", "c"])`

`T = elem_type(algebra)`

`inverses = Union{T, Tuple{T, T}}[a, b, c]`

`relators = Union{T, Tuple{T, Weight}}[(a*b)^3, (a*c)^2, (b*c)^3]`

`submodul_generators = NTuple{1, T}[(a*b*c - 1,)]`


Output:


a: \n
 0  1  0  0  0  0\n
 1  0  0  0  0  0\n
 0  0  0  1  0  0\n
 0  0  1  0  0  0\n
 0  0  0  0  0  1\n
 0  0  0  0  1  0\n

b: \n
 0  0  0  0  1  0\n
 0  0  1  0  0  0\n
 0  1  0  0  0  0\n
 0  0  0  0  0  1\n
 1  0  0  0  0  0\n
 0  0  0  1  0  0\n

c:\n
 0  0  1  0  0  0\n
 0  0  0  1  0  0\n
 1  0  0  0  0  0\n
 0  1  0  0  0  0\n
 0  0  0  0  0  1\n
 0  0  0  0  1  0\n

as sparse matrices


dimension = 6 

ntotal = 7

"""
function vector_enumeration(algebra::Any, inverses::AbstractArray, relators::AbstractArray, rank::Int=1, submodule_generators::AbstractArray=[]; abelian::Bool=false, max_weight::Int=100)
    @debug "reading input for vector enumeration over $algebra"

    if !(typeof(algebra) <: FreeAssAlgebra)
        error("`algebra` needs to be of type $FreeAssAlgebra")
    end

    gens = generators(algebra)

    A = elem_type(algebra)

    @debug "read in inverses"
    # setup inverses
    # inverses are treated specially for increased performance
    # this makes the relations xx⁻¹ = 1 and x⁻¹x = 1 redundant
    for x in inverses
        if isa(x, Tuple)
            length(x) != 2 && error("Problem with input $x: elements in `inverses` need to be either `(a, b)` with `a,b ∈ $algebra` or `x` with `x ∈ $algebra`")
            # pair of inverse generators
            x₁, x₂ = x
            (typeof(x₁) != A || typeof(x₂) != A) && error("Problem with input $x: $x in `inverses` need to be of type $A")
            !isgenerator(x₁) && error("Problem with input $x: `$x = (x₁, x₂)`, `x₁` must be a generating variable of $algebra")
            !isgenerator(x₂) && error("Problem with input $x: `$x = (x₁, x₂)`, `x₂` must be a generating variable of $algebra")
            g₁ = gens[x₁.exps[1][1]]
            g₂ = gens[x₂.exps[1][1]]

            !isnothing(g₁.inverse) && (g₁.inverse != g₂) && error("Problem with input $x: `$x₁` has already an inverse that is not `$x₂`") # there can only be one inverse to x₁
            !isnothing(g₂.inverse) && (g₂.inverse != g₁) && error("Problem with input $x: `$x₂` has already an inverse that is not `$x₁`") # there can only be one inverse to x₂
            g₁.inverse = g₂
            g₂.inverse = g₁

        elseif typeof(x) != A
            error("Problem with input $x: elements in `inverses` need to be either `(a, b)` with `a,b ∈ $algebra` or `x` with `x ∈ $algebra`")
        else
            !isgenerator(x) && error("Problem with input $x: $x must be a generating variable of $algebra")
            # involutionary / self-inverse generator
            g = gens[x.exps[1][1]]

            !isnothing(g.inverse) && (g.inverse != g) && error("Problem with input $x: `$x` has already an inverse that is not `$x`") # there can only be one inverse to x
            g.inverse = g
        end
    end

    @debug "read in relations"
    # setup relations
    rels = Vector{Relation}()

    # extra relations for abelian algebras
    if abelian
        # add relation g₁g₂ = g₂g₁ ⟺ g₁g₂ - g₂g₁ + 1 = 1 for every pair of distinct generators (g₁, g₂) apart from pairs of inverse elements
        for i = 1:length(gens)
            g₁ = gens[i]
            x₁ = element(g₁)

            for j = (i+1):length(gens)
                g₂ = gens[j]
                x₂ = element(g₂)

                if g₁.inverse != g₂
                    push!(rels, Relation(x₁ * x₂ - x₂ * x₁, gens; weight=Weight(2), fixing=false))
                end
            end
        end
    end

    # add actual relations
    for r in relators
        # weighted relators
        if isa(r, Tuple)
            length(r) != 2 && error("Problem with input $r: elements in `relators` need to be either `(a, w)` with `a ∈ $algebra` and `w` of type `Weight` or `r` with `r ∈ $algebra`")

            a, w = r
            typeof(a) != A && error("Problem with input $r: weighted relators need to be of form `(a, w)` with `a ∈ $algebra`")
            typeof(w) != Weight && error("Problem with input $r: weighted relators need to be of form `(a, w)` with `w` of type `Weight`")

            push!(rels, Relation(a, gens; weight=w))

        elseif typeof(r) != A
            error("Problem with input $r: elements in `relators` need to be either `(a, w)` with `a ∈ $algebra` and `w` of type `Weight` or `r` with `r ∈ $algebra`")

            # unweighted relators
        else
            push!(rels, Relation(r, gens))
        end
    end

    # add mandatory relations to fill the table
    for g in gens
        push!(rels, DefineRelation(g, Weight(3, 6)))
    end

    @debug "read in submodul_generators"
    # setup submodule generators
    for s in submodule_generators
        !isa(s, Tuple) && error("Problem with input $s: elements in `submodul_generators` need to be `Tuples`")
        !(length(s) == rank) && error("Problem with input $s: $s needs to be of length $rank")
        for a in s
            typeof(a) != A && error("Problem with input $s: $a needs to be in $algebra")
        end
    end
    # note: submodule generators are always interpreted as algebra words, as the overhead for missclassified group words is neglectable
    subgens = NTuple{rank,AlgebraWord}[map(s -> AlgebraWord(s, gens), tuple) for tuple in submodule_generators]


    result = run_ve(algebra, gens, rels, subgens; max_weight)

    return result
end


function run_ve(A::FreeAssAlgebra, generators::Vector{Generator}, relations::Vector{<:Relation}, subgens::Vector{NTuple{n,AlgebraWord}}; max_weight::Union{Int,Missing}=100, track_image::Bool=false) where {n}
    R = base_ring(A)
    T = elem_type(R)
    weight = 1  # start weight

    # setup the table
    @debug "creating the start basis with $n elements"
    basis = Basis(R)
    start_basis = [create!(basis, weight, length(generators), track_image ? (i, one(A)) : nothing) for i = 1:n]

    # process submodule generators at weight 1
    @debug "processing submodule generators"
    for subgen in subgens
        process_subgen!(subgen, start_basis, weight, generators)
    end

    # process relators at every weight from 2 to maximum
    weight = 2
    while isless(weight, max_weight)
        @debug "processing relations with weight $weight"

        if isclosed(basis)
            @debug "early closing at weight $weight"
            # The action on the basis elements is total, hence no new basis elements will be created by any relation.
            # Now process all relations one last time without any weight restriction to obtain equations and thus delete redundant basis elements.
            weight = missing  # note: missing is treated like ∞
        end

        #! outcommented part for lookahead !
        #=
        if activate_lookahead(basis)
            @debug "go into lookahead-mode"
            for relation in relations
                process_relation_lookahead!(basis, relation, weight; generators)
            end
            if closed(basis)
                @debug "early closing at weight $weight"
                # The action on the basis elements is total, hence no new basis elements will be created by any relation.
                # Now process all relations one last time without any weight restriction to obtain equations and thus delete redundant basis elements.
                weight = missing  # note: missing is treated like ∞
                continue
            end
        end
        =#

        for relation in relations
            process_relation!(basis, relation, weight, generators)
        end

        weight += 1
    end

    # the table should be closed by now unless the weight limit has been exceeded 
    isclosed(basis) || error("vector enumeration did not finish at maximal weight $max_weight")

    dim = nalive(basis)
    @debug "finished with dimension $dim"

    ## generate output
    # compute images of basis elements
    images = nothing
    if track_image
        @debug "computing images for output"
        images = [sparsevec([b.image[1]], [b.image[2]], n) for b in basis]
    end

    # compute preimages of standard algebra generators
    @debug "computing preimages for output"
    preimages = [SparseVector(vector(b)) for b in start_basis]

    # compute sparse matrix representation
    @debug "creating matrix representation for output"
    indices = Dict(b => i for (i, b) in enumerate(basis))  # map undeleted basis elements to their index
    matrices = Vector{SparseMatrixCSC{T,Int}}()
    for x in generators
        # generate sparse matrix
        I = Vector{Int}()
        J = Vector{Int}()
        V = Vector{T}()

        for b in basis
            @assert !isnothing(b[x])
            for t in b[x]
                c = elem(t)
                λ = coeff(t)
                push!(I, indices[b])  # row index
                push!(J, indices[c])  # column index
                push!(V, λ)  # entry
            end
        end

        push!(matrices, sparse(I, J, V, dim, dim))
    end

    @info "defined $(ncreated(basis)) basis elements to create a $dim-dimensional matrix representation"

    return (map(element, generators), matrices), (dim, ncreated(basis)), (images, preimages)  # (algebra generators, matrix generators), (dimension, total created basis elements), (images of basis elements, preimages of standard algebra generators)
end


"""
    base_qa(A::FreeAssAlgebra{FE}, R::Vector{AE}=AE[]) where {FE, AE <: FreeAssAlgElem{FE}}

Return a list of monomials of `A` forming a basis of the quotient algebra `A/⟨ARA⟩` 
of the free associative algebra `A` and the 2-sided ideal generated by the subset `R` of `A`.
"""
base_qa(A::FreeAssAlgebra{FE}, R::Vector{AE}=AE[]) where {FE,AE<:FreeAssAlgElem{FE}} = base_qm(A, 1, R, Vector{AE}[])


"""
    base_qm(A::FreeAssAlgebra{FE}, M::AbstractAlgebra.Module{AE}, R::Vector{AE}=AE[], W::Vector{ME}=ME[]) where {FE, AE <: FreeAssAlgElem{FE}, ME <: ModuleElem{AE}}

Return a list of monomials of the free module `M = Aˢ` forming a basis of the quotient module `Pˢ/⟨πˢ(W)Pˢ⟩` 
for the quotient algebra `P := A/⟨ARA⟩` of the free associative algebra `A` 
and the 2-sided ideal generated by the subset `R` of `A`, 
and the image `π(W)` of the subset `W` of `M` under the morphism `πˢ : Aˢ → Pˢ`. 

# Examples

```jldoctest; setup = :(using AbstractAlgebra, VectorEnumeration)
julia> A, (a, b, c) = free_associative_algebra(GF(7), ["a", "b", "c"])
(Free associative algebra on 3 indeterminates over finite field F_7, AbstractAlgebra.Generic.FreeAssAlgElem{AbstractAlgebra.GFElem{Int64}}[a, b, c])

julia> M = free_module(A, 1)
Free module of rank 1 over free associative algebra on 3 indeterminates over finite field F_7

julia> R = [a^2 - 1, b^2 - 1, c^2 - 1, (a*b)^3 - 1, (a*c)^2 - 1, (b*c)^3 - 1]
6-element Vector{AbstractAlgebra.Generic.FreeAssAlgElem{AbstractAlgebra.GFElem{Int64}}}:
 a^2 + 6
 b^2 + 6
 c^2 + 6
 a*b*a*b*a*b + 6
 a*c*a*c + 6
 b*c*b*c*b*c + 6

julia> W = [M([a*b*c - 1])]
1-element Vector{AbstractAlgebra.Generic.FreeModuleElem{AbstractAlgebra.Generic.FreeAssAlgElem{AbstractAlgebra.GFElem{Int64}}}}:
 (a*b*c + 6)

julia> base_qm(A, M, R, W)
6-element Vector{AbstractAlgebra.Generic.FreeModuleElem{AbstractAlgebra.Generic.FreeAssAlgElem{AbstractAlgebra.GFElem{Int64}}}}:
 (1)
 (a)
 (a*b)
 (a*b*a)
 (b)
 (a*b*a*b)
```
"""
function base_qm(A::FreeAssAlgebra{FE}, M::AbstractAlgebra.Module{AE}, R::Vector{AE}=AE[], W::Vector{ME}=ME[]) where {FE,AE<:FreeAssAlgElem{FE},ME<:ModuleElem{AE}}
    # input validation
    base_ring(M) == A || error("The free module $M is no module over $A.")

    for w in W
        parent(w) == M || error("The submodule generator $w is no element of $M.")
    end

    # return base as elements of M
    return map(M, base_qm(A, rank(M), R, map(W) do w
        [w.v...]
    end))
end


"""
    base_qm(A::FreeAssAlgebra{FE}, s::Int=1, R::Vector{AE}=AE[], W::Vector{Vector{AE}}=Vector{AE}[]) where {FE, AE <: FreeAssAlgElem{FE}}

Return a list of mononomials of `Aˢ` forming a basis of the quotient module `Pˢ/⟨πˢ(W)Pˢ⟩` 
for the quotient algebra `P := A/⟨ARA⟩` of the free associative algebra `A` 
and the 2-sided ideal generated by the subset `R` of `A`, 
and the image `π(W)` of the subset `W` of the free module `Aˢ` of rank `s` over `A`
under the morphism `πˢ : Aˢ → Pˢ`. 
"""
function base_qm(A::FreeAssAlgebra{FE}, s::Int=1, R::Vector{AE}=AE[], W::Vector{Vector{AE}}=Vector{AE}[]) where {FE,AE<:FreeAssAlgElem{FE}}
    basis = ve(A, s, R, W; track_image=true)

    @debug "collecting basis"
    # return what is stored as the image 
    z = zero(A)
    B = [[b.image[1] == i ? b.image[2] : z for i = 1:s] for b in basis]

    return B
end


"""
    dimension_qa(A::FreeAssAlgebra{FE}, R::Vector{AE}=AE[]) where {FE, AE <: FreeAssAlgElem{FE}}

Return the dimension of the quotient algebra `A/⟨ARA⟩` 
of the free associative algebra `A` and the 2-sided ideal generated by the subset `R` of `A`.
"""
dimension_qa(A::FreeAssAlgebra{FE}, R::Vector{AE}=AE[]) where {FE,AE<:FreeAssAlgElem{FE}} = dimension_qm(A, 1, R, Vector{AE}[])


"""
    dimension_qm(A::FreeAssAlgebra{FE}, M::AbstractAlgebra.Module{AE}, R::Vector{AE}=AE[], W::Vector{ME}=ME[]) where {FE, AE <: FreeAssAlgElem{FE}, ME <: ModuleElem{AE}}

Return the dimension of the quotient module `Pˢ/⟨πˢ(W)Pˢ⟩` 
for the quotient algebra `P := A/⟨ARA⟩` of the free associative algebra `A` 
and the 2-sided ideal generated by the subset `R` of `A`, 
and the image `π(W)` of the subset `W` of the free module `M = Aˢ` under the morphism `πˢ : Aˢ → Pˢ`. 

# Examples

```jldoctest; setup = :(using AbstractAlgebra, VectorEnumeration)
julia> A, (a, b, c) = free_associative_algebra(GF(7), ["a", "b", "c"])
(Free associative algebra on 3 indeterminates over finite field F_7, AbstractAlgebra.Generic.FreeAssAlgElem{AbstractAlgebra.GFElem{Int64}}[a, b, c])

julia> M = free_module(A, 1)
Free module of rank 1 over free associative algebra on 3 indeterminates over finite field F_7

julia> R = [a^2 - 1, b^2 - 1, c^2 - 1, (a*b)^3 - 1, (a*c)^2 - 1, (b*c)^3 - 1]
6-element Vector{AbstractAlgebra.Generic.FreeAssAlgElem{AbstractAlgebra.GFElem{Int64}}}:
 a^2 + 6
 b^2 + 6
 c^2 + 6
 a*b*a*b*a*b + 6
 a*c*a*c + 6
 b*c*b*c*b*c + 6

julia> W = [M([a*b*c - 1])]
1-element Vector{AbstractAlgebra.Generic.FreeModuleElem{AbstractAlgebra.Generic.FreeAssAlgElem{AbstractAlgebra.GFElem{Int64}}}}:
 (a*b*c + 6)

julia> dimension_qm(A, M, R, W)
6
```
"""
function dimension_qm(A::FreeAssAlgebra{FE}, M::AbstractAlgebra.Module{AE}, R::Vector{AE}=AE[], W::Vector{ME}=ME[]) where {FE,AE<:FreeAssAlgElem{FE},ME<:ModuleElem{AE}}
    # input validation
    base_ring(M) == A || error("The free module $M is no module over $A.")

    for w in W
        parent(w) == M || error("The submodule generator $w is no element of $M.")
    end

    return dimension_qm(A, rank(M), R, map(W) do w
        [w.v...]
    end)
end


"""
    dimension_qm(A::FreeAssAlgebra{FE}, s::Int=1, R::Vector{AE}=AE[], W::Vector{Vector{AE}}=Vector{AE}[]) where {FE, AE <: FreeAssAlgElem{FE}}

Return the dimension of the quotient module `Pˢ/⟨πˢ(W)Pˢ⟩` 
for the quotient algebra `P := A/⟨ARA⟩` of the free associative algebra `A` 
and the 2-sided ideal generated by the subset `R` of `A`, 
and the image `π(W)` of the subset `W` of the free module `Aˢ` of rank `s` over `A`
under the morphism `πˢ : Aˢ → Pˢ`. 
"""
function dimension_qm(A::FreeAssAlgebra{FE}, s::Int=1, R::Vector{AE}=AE[], W::Vector{Vector{AE}}=Vector{AE}[]) where {FE,AE<:FreeAssAlgElem{FE}}
    basis = ve(A, s, R, W)
    dim = nalive(basis)

    return dim
end


"""
    matrices_qa(T::Type, A::FreeAssAlgebra{FE}, R::Vector{AE}=AE[]) where {FE, AE <: FreeAssAlgElem{FE}}

Return matrix representations as matrices of type `T` for the generators of the free associative algebra `A` for right multiplication in the quotient algebra `A/⟨ARA⟩` 
of `A` and the 2-sided ideal generated by the subset `R` of `A`.
"""
matrices_qa(T::Type, A::FreeAssAlgebra{FE}, R::Vector{AE}=AE[]) where {FE,AE<:FreeAssAlgElem{FE}} = matrices_qm(T, A, 1, R, Vector{AE}[])


"""
    matrices_qm(T::Type, A::FreeAssAlgebra{FE}, M::AbstractAlgebra.Module{AE}, R::Vector{AE}=AE[], W::Vector{ME}=ME[]) where {FE, AE <: FreeAssAlgElem{FE}, ME <: ModuleElem{AE}}

Return matrix representations as matrices of type `T` for the generators of the free associative algebra `A` for right multiplication in the quotient module `Pˢ/⟨πˢ(W)Pˢ⟩` 
for the quotient algebra `P := A/⟨ARA⟩` of `A` and the 2-sided ideal generated by the subset `R` of `A`, 
and the image `π(W)` of the subset `W` of the free module `M = Aˢ` under the morphism `πˢ : Aˢ → Pˢ`. 

# Examples

```jldoctest; setup = :(using AbstractAlgebra, VectorEnumeration)
julia> A, (a, b, c) = free_associative_algebra(GF(7), ["a", "b", "c"])
(Free associative algebra on 3 indeterminates over finite field F_7, AbstractAlgebra.Generic.FreeAssAlgElem{AbstractAlgebra.GFElem{Int64}}[a, b, c])

julia> M = free_module(A, 1)
Free module of rank 1 over free associative algebra on 3 indeterminates over finite field F_7

julia> R = [a^2 - 1, b^2 - 1, c^2 - 1, (a*b)^3 - 1, (a*c)^2 - 1, (b*c)^3 - 1]
6-element Vector{AbstractAlgebra.Generic.FreeAssAlgElem{AbstractAlgebra.GFElem{Int64}}}:
 a^2 + 6
 b^2 + 6
 c^2 + 6
 a*b*a*b*a*b + 6
 a*c*a*c + 6
 b*c*b*c*b*c + 6

julia> W = [M([a*b*c - 1])]
1-element Vector{AbstractAlgebra.Generic.FreeModuleElem{AbstractAlgebra.Generic.FreeAssAlgElem{AbstractAlgebra.GFElem{Int64}}}}:
 (a*b*c + 6)

julia> Xᵃ, Xᵇ, Xᶜ = matrices_qm(Matrix, A, M, R, W)
3-element Vector{Matrix{AbstractAlgebra.GFElem{Int64}}}:
 [0 1 … 0 0; 1 0 … 0 0; … ; 0 0 … 0 1; 0 0 … 1 0]
 [0 0 … 1 0; 0 0 … 0 0; … ; 1 0 … 0 0; 0 0 … 0 0]
 [0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 1; 0 0 … 1 0]

julia> Xᵃ
6×6 Matrix{AbstractAlgebra.GFElem{Int64}}:
 0  1  0  0  0  0
 1  0  0  0  0  0
 0  0  0  1  0  0
 0  0  1  0  0  0
 0  0  0  0  0  1
 0  0  0  0  1  0

julia> Xᵇ
6×6 Matrix{AbstractAlgebra.GFElem{Int64}}:
 0  0  0  0  1  0
 0  0  1  0  0  0
 0  1  0  0  0  0
 0  0  0  0  0  1
 1  0  0  0  0  0
 0  0  0  1  0  0

julia> Xᶜ
6×6 Matrix{AbstractAlgebra.GFElem{Int64}}:
 0  0  1  0  0  0
 0  0  0  1  0  0
 1  0  0  0  0  0
 0  1  0  0  0  0
 0  0  0  0  0  1
 0  0  0  0  1  0
```
"""
function matrices_qm(T::Type, A::FreeAssAlgebra{FE}, M::AbstractAlgebra.Module{AE}, R::Vector{AE}=AE[], W::Vector{ME}=ME[]) where {FE,AE<:FreeAssAlgElem{FE},ME<:ModuleElem{AE}}
    # input validation
    base_ring(M) == A || error("The free module $M is no module over $A.")

    for w in W
        parent(w) == M || error("The submodule generator $w is no element of $M.")
    end

    return matrices_qm(T, A, rank(M), R, map(W) do w
        [w.v...]
    end)
end

"""
    matrices_qm(T::Type, A::FreeAssAlgebra{FE}, s::Int=1, R::Vector{AE}=AE[], W::Vector{Vector{AE}}=Vector{AE}[]) where {FE, AE <: FreeAssAlgElem{FE}}

Return matrix representations as matrices of type `T` for the generators of the free associative algebra `A` for right multiplication in the quotient module `Pˢ/⟨πˢ(W)Pˢ⟩` 
for the quotient algebra `P := A/⟨ARA⟩` of `A` and the 2-sided ideal generated by the subset `R` of `A`, 
and the image `π(W)` of the subset `W` of the free module `Aˢ` of rank `s` over `A` under the morphism `πˢ : Aˢ → Pˢ`. 
"""
function matrices_qm(::Type{SparseMatrixCSC}, A::FreeAssAlgebra{FE}, s::Int=1, R::Vector{AE}=AE[], W::Vector{Vector{AE}}=Vector{AE}[]) where {FE,AE<:FreeAssAlgElem{FE}}
    basis = ve(A, s, R, W)
    dim = nalive(basis)

    @debug "creating matrix representations"
    indices = Dict(b => i for (i, b) in enumerate(basis))  # map undeleted basis elements to their index
    matrices = Vector{SparseMatrixCSC{FE,Int}}()
    for x = 1:ngens(A)
        # generate sparse matrix
        I = Vector{Int}()
        J = Vector{Int}()
        V = Vector{FE}()

        for b in basis
            @assert !isnothing(b[x])
            for t in b[x]
                c = elem(t)
                λ = coeff(t)
                push!(I, indices[b])  # row index
                push!(J, indices[c])  # column index
                push!(V, λ)  # entry
            end
        end

        push!(matrices, sparse(I, J, V, dim, dim))
    end

    return matrices
end

function matrices_qm(::Type{Matrix}, A::FreeAssAlgebra{FE}, s::Int=1, R::Vector{AE}=AE[], W::Vector{Vector{AE}}=Vector{AE}[]) where {FE,AE<:FreeAssAlgElem{FE}}
    basis = ve(A, s, R, W)
    dim = nalive(basis)

    @debug "creating matrix representations"
    indices = Dict(b => i for (i, b) in enumerate(basis))  # map undeleted basis elements to their index
    matrices = Vector{Matrix{FE}}()
    z = zero(base_ring(A))
    for x = 1:ngens(A)
        # generate dense matrix
        M = fill(z, (dim, dim))

        for b in basis

            @assert !isnothing(b[x])
            for t in b[x]
                c = elem(t)
                λ = coeff(t)
                M[indices[b], indices[c]] = λ
            end
        end

        push!(matrices, M)
    end

    return matrices
end


function ve(A::FreeAssAlgebra{FE}, s::Int, R::Vector{AE}, W::Vector{Vector{AE}}; max_weight::Union{Int,Missing}=100, track_image::Bool=false) where {FE,AE<:FreeAssAlgElem{FE}}
    # rank must be positive
    s > 0 || error("The module rank $s must be positive.")

    # setup generators
    gens = [Generator(i, gen(A, i)) for i in 1:nvars(A)]

    # setup relations
    R = Set(R) # drop duplicate relators
    relations = Relation[]

    for r in R
        # validate relators are elements of the algebra
        parent(r) == A || error("The relator $r is no element of $A.")

        # relations for zero relators are trivially satisfied 
        iszero(r) && continue

        # extract relations for inverse generators, i.e. relators of type x*y - 1 for generators x, y for which also the relator y*x - 1 is present
        if length(r) == 2
            m₁ = exponent_word(r, 1)
            m₂ = exponent_word(r, 2)

            # monomials are sorted by degree in free associative algebra elements, so deg(m₁) ≥ deg(m₂)
            if length(m₁) == 2 && length(m₂) == 0
                # relation x*y - 1
                x, y = m₁

                # TODO: what to do, if there are multiple inverses, e.g. if there are duplicate generators?
                # TODO: add warning message 
                if !is_invertible(gens[x]) && !is_invertible(gens[y]) && (x == y || (gen(A, y) * gen(A, x) - 1) in R)
                    @debug "found inverse $(gen(A, y)) to generator $(gen(A, x))"
                    gens[x].inverse = gens[y]
                    gens[y].inverse = gens[x]
                    continue  # exclude the relation obtained from this relator
                end
            end
        end

        # add relation
        # TODO: inline Relation and work with annihilating relators to begin with
        push!(relations, Relation(r + 1, gens))
    end

    # add virtual relations x = x for every generator x to guarantee the action is total on termination
    for x in gens
        push!(relations, DefineRelation(x, Weight(3, 6)))
    end

    # setup submodule generators
    # note: submodule generators are always interpreted as algebra words, as the overhead for missclassified group words is neglectable
    # TODO: maybe it would be just better not to use tuples here
    subgens = NTuple{s,AlgebraWord}[]
    for w in W
        # validate all submodule generators have the correct rank and their entries are elements of the algebra
        length(w) == s || error("The submodule generator $w is not of rank $s.")
        for a in w
            parent(a) == A || error("The entry $a of the submodule generator $w is no element of $A.")
        end

        push!(subgens, NTuple{s,AlgebraWord}(AlgebraWord(a, gens) for a in w))
    end

    ### run vector enumeration
    # setup the start basis 
    weight = 1  # start weight

    @debug "creating the start basis with $s elements"
    basis = Basis(base_ring(A))
    start_basis = [create!(basis, weight, length(gens), track_image ? (i, one(A)) : nothing) for i = 1:s]

    # process submodule generators at weight 1
    @debug "processing submodule generators"
    for subgen in subgens
        process_subgen!(subgen, start_basis, weight, gens)
    end

    # process relators at every weight from 2 to maximum
    weight = 2
    while isless(weight, max_weight)
        @debug "processing relations with weight $weight"

        if isclosed(basis)
            @debug "early closing at weight $weight"
            # The action on the basis elements is total, hence no new basis elements will be created by any relation.
            # Now process all relations one last time without any weight restriction to obtain equations and thus delete redundant basis elements.
            weight = missing  # note: missing is treated like ∞
        end

        #! outcommented part for lookahead !
        #=
        if activate_lookahead(basis)
            @debug "go into lookahead-mode"
            for relation in relations
                process_relation_lookahead!(basis, relation, weight; generators)
            end
            if closed(basis)
                @debug "early closing at weight $weight"
                # The action on the basis elements is total, hence no new basis elements will be created by any relation.
                # Now process all relations one last time without any weight restriction to obtain equations and thus delete redundant basis elements.
                weight = missing  # note: missing is treated like ∞
                continue
            end
        end
        =#

        for relation in relations
            process_relation!(basis, relation, weight, gens)
        end

        weight += 1
    end

    # the table should be closed by now unless the weight limit has been exceeded 
    isclosed(basis) || error("vector enumeration did not finish at maximal weight $max_weight")

    dim = nalive(basis)
    @debug "finished with dimension $dim after defining $(ncreated(basis)) basis elements"

    # return the basis state from which all output can be obtained (except for the start basis (the preimages), which is not exported right now)
    return basis
end

end # module VectorEnumeration
