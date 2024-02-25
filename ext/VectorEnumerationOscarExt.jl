module VectorEnumerationOscarExt

using VectorEnumeration, Oscar

import VectorEnumeration: base_qa, base_qm, dimension_qa, dimension_qm, matrices_qa, matrices_qm

#TODO: On further development of submodules in AbstractAlgebra change S to Submodule of M

"""
    base_qa(A::FreeAssAlgebra{FE}, I::Oscar.FreeAssAlgIdeal{AE}) where {FE, AE <: FreeAssAlgElem{FE}}

Return a subset of `A` forming a k basis of the quotient algebra `A/I` 
of the free associative algebra `A` and the 2-sided ideal `I`.
"""
base_qa(A::FreeAssAlgebra{FE}, I::Oscar.FreeAssAlgIdeal{AE}) where {FE, AE <: FreeAssAlgElem{FE}} = base_ring(I) == A ? base_qa(T, A, gens(I)) : error("Ideal $I needs to be in $A")

"""
    base_qm(A::FreeAssAlgebra{FE}, M::AbstractAlgebra.Module{AE}, I::Oscar.FreeAssAlgIdeal{AE}, W::Vector{ME}=ME[]) where {FE, AE <: FreeAssAlgElem{FE}, ME <: ModuleElem{AE}}

Return a subset of the free module `M = Aˢ` forming a basis of the quotient module `Pˢ/⟨πˢ(W)Pˢ⟩` 
for the quotient algebra `P := A/I` of the free associative algebra `A` 
and the 2-sided ideal `I`, 
and the image `π(W)` of the subset `W` of `M` under the morphism `πˢ : Aˢ → Pˢ`. 
"""
function base_qm(A::FreeAssAlgebra{FE}, M::AbstractAlgebra.Module{AE}, I::Oscar.FreeAssAlgIdeal{AE}, W::Vector{ME}=ME[]) where {FE, AE <: FreeAssAlgElem{FE}, ME <: ModuleElem{AE}}
    base_ring(I) == A || error("Ideal $I needs to be in $A")

    return base_qm(A, M, gens(I), W)
end

"""
    dimension_qa(A::FreeAssAlgebra{FE}, I::Oscar.FreeAssAlgIdeal{AE}) where {FE, AE <: FreeAssAlgElem{FE}}

Return the dimension of the quotient algebra `A/I` 
of the free associative algebra `A` and the 2-sided ideal `I`.
"""
dimension_qa(A::FreeAssAlgebra{FE}, I::Oscar.FreeAssAlgIdeal{AE}) where {FE, AE <: FreeAssAlgElem{FE}} = base_ring(I) == A ? dimension_qa(T, A, gens(I)) : error("Ideal $I needs to be in $A")

"""
    dimension_qm(A::FreeAssAlgebra{FE}, M::AbstractAlgebra.Module{AE}, I::Oscar.FreeAssAlgIdeal{AE}, W::Vector{ME}=ME[]) where {FE, AE <: FreeAssAlgElem{FE}, ME <: ModuleElem{AE}}

Return the dimension of the quotient module `Pˢ/⟨πˢ(W)Pˢ⟩` 
for the quotient algebra `P := A/I` of the free associative algebra `A` 
and the 2-sided ideal `I`, 
and the image `π(W)` of the subset `W` of the free module `M = Aˢ` under the morphism `πˢ : Aˢ → Pˢ`. 
"""
function dimension_qm(A::FreeAssAlgebra{FE}, M::AbstractAlgebra.Module{AE}, I::Oscar.FreeAssAlgIdeal{AE}, W::Vector{ME}=ME[]) where {FE, AE <: FreeAssAlgElem{FE}, ME <: ModuleElem{AE}}
    base_ring(I) == A || error("Ideal $I needs to be in $A")

    return dimension_qm(A, M, gens(I), W)
end

"""
    matrices_qa(T::Type, A::FreeAssAlgebra{FE}, I::Oscar.FreeAssAlgIdeal{AE}) where {FE, AE <: FreeAssAlgElem{FE}}

Return matrix representations as matrices of type `T` for the generators of the free associative algebra `A` for right multiplication in the quotient algebra `A/I` 
of `A` and the 2-sided ideal `I`.
"""
matrices_qa(T::Type, A::FreeAssAlgebra{FE}, I::Oscar.FreeAssAlgIdeal{AE}) where {FE, AE <: FreeAssAlgElem{FE}} = base_ring(I) == A ? matrices_qa(T, A, gens(I)) : error("Ideal $I needs to be in $A")

"""
    matrices_qm(T::Type, A::FreeAssAlgebra{FE}, M::AbstractAlgebra.Module{AE}, I::Oscar.FreeAssAlgIdeal{AE}, W::Vector{ME}=ME[]) where {FE, AE <: FreeAssAlgElem{FE}, ME <: ModuleElem{AE}}

Return matrix representations as matrices of type `T` for the generators of the free associative algebra `A` for right multiplication in the quotient module `Pˢ/⟨πˢ(W)Pˢ⟩` 
for the quotient algebra `P := A/I` of `A` and the 2-sided ideal `I`, 
and the image `π(W)` of the subset `W` of the free module `M = Aˢ` under the morphism `πˢ : Aˢ → Pˢ`. 
"""
function matrices_qm(T::Type, A::FreeAssAlgebra{FE}, M::AbstractAlgebra.Module{AE}, I::Oscar.FreeAssAlgIdeal{AE}, W::Vector{ME}=ME[]) where {FE, AE <: FreeAssAlgElem{FE}, ME <: ModuleElem{AE}}
    base_ring(I) == A || error("Ideal $I needs to be in $A")

    return matrices_qm(T, A, M, gens(I), W)
end


function setupgrouprelation(f::Field, FP::Oscar.FPGroup, U::Oscar.FPGroup)
    n = ngens(FP)

    relations = relators(FP)
    
    involutions = []

    relations = map(x -> map_word(x, [FP[i] for i in 1:n]), relations) #TODO find better way

    for i in 1:n
        if FP[i]^2 in relations 
            push!(involutions, i)
        end
    end
    
    A, G = free_associative_algebra(f, 2*n - length(involutions))
    
    genimg = G[1:n]
    genimginv = elem_type(A)[]

    i₀ = n + 1
    for i in 1:n
        if i in involutions
            push!(genimginv, G[i])
        else
            push!(genimginv, G[i₀])
            i₀ += 1
        end
    end

    
    relations = map(x -> map_word(x, genimg; genimgs_inv = genimginv) - one(A), relations)
    
    for i in 1:n
        if !(i in involutions)
            push!(relations, genimg[i]*genimginv[i] - one(A))
        end
    end
    
    submodule_generators = map(x -> [map_word(x, genimg; genimgs_inv = genimginv) - one(A)], gens(U))

    return A, relations, submodule_generators
end

"""
    base_qm(f::Field, FP::Oscar.FPGroup, U::Oscar.FPGroup=trivial_subgroup(FP))

Return a list of monomials forming a f - basis of f(FP/U)

# Examples

```jldoctest; setup = :(using Oscar, VectorEnumeration)
julia> F = free_group(3)
<free group on the generators [ f1, f2, f3 ]>

julia> a, b, c = gens(F)
3-element Vector{FPGroupElem}:
 f1
 f2
 f3

julia> R = [a^2, b^2, c^2, (a*b)^3, (a*c)^2, (b*c)^3]
6-element Vector{FPGroupElem}:
 f1^2
 f2^2
 f3^2
 (f1*f2)^3
 (f1*f3)^2
 (f2*f3)^3

julia> N = normal_closure(F, sub(F, R)[1])[1]
Group(<free, no generators known>)

julia> FP = quo(FPGroup, F, N)[1]
<fp group of size 24 on the generators [ F1, F2, F3 ]>

julia> U = sub(FP, [FP[1]*FP[2]*FP[3]])[1]
Group([ F1*F2*F3 ])

julia> base_qm(QQ, FP, U)
6-element Vector{Vector{AbstractAlgebra.Generic.FreeAssAlgElem{QQFieldElem}}}:
 [1]
 [x1]
 [x1*x2]
 [x1*x2*x1]
 [x2]
 [x1*x2*x1*x2]
```
"""
function base_qm(f::Field, FP::Oscar.FPGroup, U::Oscar.FPGroup=trivial_subgroup(FP))

    is_full_fp_group(FP) || error("$FP needs to be a finitely presented group")

    A, relations, submodule_generators = setupgrouprelation(f, FP, U)

    return base_qm(A, 1, relations, submodule_generators)
end


"""
    dimension_qm(f::Field, FP::Oscar.FPGroup, U::Oscar.FPGroup=trivial_subgroup(FP))

Return the f - dimension of f(FP/U)

# Examples

```jldoctest; setup = :(using Oscar, VectorEnumeration)
julia> F = free_group(3)
<free group on the generators [ f1, f2, f3 ]>

julia> a, b, c = gens(F)
3-element Vector{FPGroupElem}:
 f1
 f2
 f3

julia> R = [a^2, b^2, c^2, (a*b)^3, (a*c)^2, (b*c)^3]
6-element Vector{FPGroupElem}:
 f1^2
 f2^2
 f3^2
 (f1*f2)^3
 (f1*f3)^2
 (f2*f3)^3

julia> N = normal_closure(F, sub(F, R)[1])[1]
Group(<free, no generators known>)

julia> FP= quo(FPGroup, F, N)[1]
<fp group of size 24 on the generators [ F1, F2, F3 ]>

julia> U = sub(FP, [FP[1]*FP[2]*FP[3]])[1]
Group([ F1*F2*F3 ])

julia> dimension_qm(QQ, FP, U)
6
```
"""
function dimension_qm(f::Field, FP::Oscar.FPGroup, U::Oscar.FPGroup=trivial_subgroup(FP))

    is_full_fp_group(FP) || error("$FP needs to be a finitely presented group")

    A, relations, submodule_generators = setupgrouprelation(f, FP, U)

    return dimension_qm(A, 1, relations, submodule_generators)
end

"""
    matrices_qm(T::Type, f::Field, FP::Oscar.FPGroup, U::Oscar.FPGroup=trivial_subgroup(FP))

Return a f - matrix - representation of f(FP/U) in matrices of type T

# Examples

```jldoctest; setup = :(using Oscar, VectorEnumeration)
julia> F = free_group(3)
<free group on the generators [ f1, f2, f3 ]>

julia> a, b, c = gens(F)
3-element Vector{FPGroupElem}:
 f1
 f2
 f3

julia> R = [a^2, b^2, c^2, (a*b)^3, (a*c)^2, (b*c)^3]
6-element Vector{FPGroupElem}:
 f1^2
 f2^2
 f3^2
 (f1*f2)^3
 (f1*f3)^2
 (f2*f3)^3

julia> N = normal_closure(F, sub(F, R)[1])[1]
Group(<free, no generators known>)

julia> FP= quo(FPGroup, F, N)[1]
<fp group of size 24 on the generators [ F1, F2, F3 ]>

julia> U = sub(FP, [FP[1]*FP[2]*FP[3]])[1]
Group([ F1*F2*F3 ])

julia> Xᵃ, Xᵇ, Xᶜ = matrices_qm(Matrix, QQ, FP, U)
3-element Vector{Matrix{QQFieldElem}}:
 [0 1 … 0 0; 1 0 … 0 0; … ; 0 0 … 0 1; 0 0 … 1 0]
 [0 0 … 1 0; 0 0 … 0 0; … ; 1 0 … 0 0; 0 0 … 0 0]
 [0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 1; 0 0 … 1 0]

julia> Xᵃ
6×6 Matrix{QQFieldElem}:
 0  1  0  0  0  0
 1  0  0  0  0  0
 0  0  0  1  0  0
 0  0  1  0  0  0
 0  0  0  0  0  1
 0  0  0  0  1  0

julia> Xᵇ
6×6 Matrix{QQFieldElem}:
 0  0  0  0  1  0
 0  0  1  0  0  0
 0  1  0  0  0  0
 0  0  0  0  0  1
 1  0  0  0  0  0
 0  0  0  1  0  0

julia> Xᶜ
6×6 Matrix{QQFieldElem}:
 0  0  1  0  0  0
 0  0  0  1  0  0
 1  0  0  0  0  0
 0  1  0  0  0  0
 0  0  0  0  0  1
 0  0  0  0  1  0
```
"""
function matrices_qm(T::Type, f::Field, FP::Oscar.FPGroup, U::Oscar.FPGroup=trivial_subgroup(FP))

    is_full_fp_group(FP) || error("$FP needs to be a finitely presented group")

    A, relations, submodule_generators = setupgrouprelation(f, FP, U)

    return matrices_qm(T, A, 1, relations, submodule_generators)    
end

end  # module VectorEnumerationOscarExt
