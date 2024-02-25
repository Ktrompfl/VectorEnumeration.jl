@testset "vector.jl" begin
    field = GF(13)
    x = Symbol("x")
    A, X = free_associative_algebra(field, 5, x)
    gen = VectorEnumeration.generators(A)
    basis = VectorEnumeration.Basis(field)
    word = VectorEnumeration.GroupWord(A(1), gen)
    B = [VectorEnumeration.create!(basis, VectorEnumeration.BasisElementData(1, word, 5, 1))]
    for i in 1:5
        word = word*gen[i]
        push!(B,  VectorEnumeration.create!(basis, VectorEnumeration.BasisElementData(i + 1, word, 5, 2)))
    end

    b = basis.last
    # B = [b₁, b₂, b₃, b₄, b₅, b₆]


    v = VectorEnumeration.zero(b) #0
    w = VectorEnumeration.SparseBasisVector(b, field(3)) #3b₆
    x = VectorEnumeration.zero(w)
    a = VectorEnumeration.zero(w)


    b = basis.first_undeleted #b₁
    VectorEnumeration.add_scaled!(x, b, field(4))
    VectorEnumeration.add_scaled!(a, b, field(4))
    b = VectorEnumeration.next(b) #b₂
    VectorEnumeration.add_scaled!(x, b, field(3)) #x = 4b₁ + 3b₂
    VectorEnumeration.add_scaled!(a, b, field(3))

    b = VectorEnumeration.next(b) #b₃

    y = x + B[2] # 4b₁ + 4b₂

   

    VectorEnumeration.delete!(basis.last, x) #delete b₆ replacing with 4b₁ + 3b₂

    z = y + w # 3b₁


    w = VectorEnumeration.SparseBasisVector(basis.last, field(3)) #3b₆

    VectorEnumeration.add_scaled!(w, a, field(1))
    VectorEnumeration.add_scaled!(w, B[5], field(7)) # 4b₁ + 3b₂ + 7b₅ + 3b₆

    VectorEnumeration.delete!(b, z) #delete b₃ replacing with 8b₁ + 7b₂
    

    @testset "Basis Tests and " begin
        @test basis.deleted == 2
        @test basis.created == 6
        @test basis.first_undeleted == B[1]
        @test basis.last_undeleted == B[5]
        @test basis.last == B[6]
        @test VectorEnumeration.base_ring(basis) === field
        l = [1,2,4,5]
        k = 1
        for b in basis
            @test !VectorEnumeration.isdeleted(b)
            @test b == B[l[k]]
            k += 1
        end
        @test (k-1) == VectorEnumeration.nalive(basis)
    end
    @testset "BasisElement Test" begin
        @test VectorEnumeration.replacement(B[3]) === z && VectorEnumeration.replacement(B[6]) === x 
        @test VectorEnumeration.next(B[2]) === B[4]
    end
    @testset "SparseVector Test" begin
        @testset "Vector properties" begin
            @testset "v" begin
                @test VectorEnumeration.basis(v) == basis
                @test isempty(v.coeff)
                @test isempty(v.elem)
                @test iszero(v)
            end
            @testset "w" begin
                @test VectorEnumeration.basis(w) == basis
                @test w.coeff == [field(3), field(7), field(3), field(4)]
                @test w.elem == [B[6], B[5], B[2],B[1]]
            end
        end

        @testset "Arithmetics" begin
            VectorEnumeration.add_scaled!(w, v, field(2))
            @test w.coeff == [field(3), field(7), field(3), field(4)]  
            @test w.elem == [B[6], B[5], B[2],B[1]]                    
            VectorEnumeration.add_scaled!(z, w, field(4)) # 2b₁ + 9b₂ + 2b₅
            @test z.coeff == [field(2), field(9), field(2)]
            @test z. elem == [B[5], B[2], B[1]]
            VectorEnumeration.add_scaled!(y, B[2], field(10)) #4b₁ + 1b₂
            @test y.coeff == [field(1), field(4)]
            @test y.elem == [B[2], B[1]]
        end
        
        
        @testset "Iteration" begin
            l = [5, 2, 1]
            c = [7, 12, 3]
            k = 1
            for (b, λ) in w  #3b₁ + 12b₂ + 7b₅ 
                @test b == B[l[k]]
                @test λ == field(c[k])
                k += 1
            end
        end
    end
end
