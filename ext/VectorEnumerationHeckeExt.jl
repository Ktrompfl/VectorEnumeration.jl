module VectorEnumerationHeckeExt

using VectorEnumeration, Hecke

import VectorEnumeration: matrices_qm

function matrices_qm(::Type{SMat}, A::FreeAssAlgebra{FE}, s::Int=1, R::Vector{AE}=AE[], W::Vector{Vector{AE}}=Vector{AE}[]) where {FE, AE <: FreeAssAlgElem{FE}}
    basis = VectorEnumeration.ve(A, s, R, W)
    dim = VectorEnumeration.nalive(basis)

    @debug "creating matrix representations"
    indices = Dict(b => i for (i, b) in enumerate(basis))  # map undeleted basis elements to their index

    matrices = Vector{SMat{FE}}()
    field = base_ring(A)
    for x = 1:ngens(A)
        # create hecke sparse matrix (row major)
        M = sparse_matrix(field, dim, dim)

        for b in basis
            @assert !isnothing(b[x])
            for t in b[x]
                c = VectorEnumeration.elem(t)
                λ = VectorEnumeration.coeff(t)
                row = M[indices[c]]
                push!(row.pos, indices[b])
                push!(row.values, λ)
                M.nnz += 1
            end
        end

        push!(matrices, M)
    end

    return matrices
end

end  # module VectorEnumerationHeckeExt