function verify(relations, matrices::Vector{SparseMatrixCSC{Tv, Ti}}) where {Tv, Ti}
    #if there are no relations to verify return
    if isempty(relations)
        return
    end

    n, m = size(matrices[1])
    @assert n == m

    # skip if dimension too large, as matrix computations get extremly expensive
    if n > 100
        return
    end

    #check relations as sparse matrices
    result = spzeros(n,n)

    #evaluate relations and check if they are fulfilled
    for r in relations
        for (c, m) in zip(r.coeffs, r.exps)
            tmp = sparse(I, n, n)
            for i in m
                tmp = tmp*matrices[i]
            end
            result = result + c*tmp
        end
        @testset "Relation $r" begin
            @test result ==  spzeros(n,n)
        end
    end
end

function verify(relations, matrices::Vector{Matrix{T}}) where {T}
    #if there are no relations to verify return
    if isempty(relations)
        return
    end

    n, m = size(matrices[1])

    # skip if dimension too large, as matrix computations get extremly expensive
    if n > 100
        return
    end

    @assert n == m

    #check relations on matrices in MatSpace 
    M = matrix_space(base_ring(parent(relations[1])), n, n)

    result = zero(M)

    #evaluate relations and check if they are fulfilled
    for r in relations
        for (c, m) in zip(r.coeffs, r.exps)
            tmp = M(Matrix(I, n, n))
            for i in m
                tmp = tmp*M(matrices[i])
            end
            result = result + c*tmp
        end
        @testset "Relation $r" begin
            @test result ==  zero(M)
        end
    end       
end
