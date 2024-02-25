import Nemo: CalciumField

@testset "algerbra.jl" begin
    @testset "Generator" begin
        x = Symbol("x")
        fields = [QQ, GF(3), CalciumField(), GF(11), GF(47)]
        for i in 1:5
            field = fields[i]
            n = rand(1:50)
            @testset "Test $i: Algebra in $n variables over $field" begin
                A ,X = free_associative_algebra(field, n, x)
                gen = VectorEnumeration.generators(A)
                @testset "isgenerator" begin
                    @testset "base_ring" begin
                        for k in 1:n
                            @test VectorEnumeration.base_ring(gen[k]) == field
                        end
                    end
                    @testset "generators" begin
                        for k in 1:n
                            @test VectorEnumeration.isgenerator(gen[k].element)
                        end    
                    end
                    @testset "non-generators" begin
                        @test !VectorEnumeration.isgenerator(gen[1].element^2)
                    end
                end
                @testset "inverse" begin
                    for k in 1:n
                        @test !VectorEnumeration.is_invertible(gen[k])
                    end
                    m = rand(1:n)
                    l = rand(1:n)
                    gen[m].inverse = gen[l]
                    gen[l].inverse = gen[m]
                    @test VectorEnumeration.is_invertible(gen[m]) && VectorEnumeration.is_invertible(gen[l])
                    @test VectorEnumeration.inverse(gen[m]) == gen[l] && VectorEnumeration.inverse(gen[l]) == gen[m]
                end
            end
        end
    end
    @testset "GroupWord" begin
        x = Symbol("x")
        fields = [QQ, GF(3), CalciumField(), GF(11), GF(47)]
        for i in 1:5
            field = fields[i]
            n = rand(1:50)
            @testset "Test $i: Algebra in $n variables over $field" begin
                A ,X = free_associative_algebra(field, n, x)
                gen = VectorEnumeration.generators(A)
                elem = A(1)
                a = A(1)
                @testset "isgroupword" begin
                    for r in 1:10
                        for k in 1:rand(1:10)
                            exp = rand(1:25)
                            elem *= gen[rand(1:n)].element^exp
                        end
                        a += elem
                        @test VectorEnumeration.isgroupword(elem)
                        @test !VectorEnumeration.isgroupword(a) || r == 1 
                    end
                end
                @testset "GroupWord()" begin 
                    for r in 1:10
                        elem = A(1)
                        for k in 1:rand(1:25)
                            exp = rand(1:10)
                            idx = rand(1:n)
                            elem = elem * (gen[idx].element)^exp
                            if k == 1 #! Vector{VectorEnumeration.Generator}() not working
                                global word = fill(gen[idx], exp)
                            else
                                append!(word, fill(gen[idx], exp))
                            end
                            
                        end
                        g = VectorEnumeration.GroupWord(elem, gen)
                        @test VectorEnumeration.base_ring(g) == field
                        @test g.element == elem
                        @test g.word == word 
                    end
                end
            end
        end
    end
    @testset "AlgebraWord" begin
        x = Symbol("x")
        fields = [QQ, GF(3), CalciumField(), GF(11), GF(47)]
        for i in 1:5
            field = fields[i]
            n = rand(1:50)
            @testset "Test $i: Algebra in $n variables over $field" begin
                A ,X = free_associative_algebra(field, n, x)
                gen = VectorEnumeration.generators(A)
                a = A(0)
                ln = rand(5:5)
                for r in 1:ln
                    elem = A(1)
                    for k in 1:rand(1:10)
                        exp = rand(1:10)
                        elem *= gen[rand(1:n)].element^exp
                    end
                    coeff = field(rand(-50:50))
                    while isequal(coeff, field(0))
                        coeff = field(rand(-50:50))
                    end
                    if r ==1
                        global coeffs = [coeff]
                        global monom = [VectorEnumeration.GroupWord(elem, gen)]
                    else
                        push!(coeffs, coeff)
                        push!(monom, VectorEnumeration.GroupWord(elem, gen))
                    end
                    a += coeff*elem
                end
                algword = VectorEnumeration.AlgebraWord(a, gen)
                @test VectorEnumeration.base_ring(algword) == field
                @test Set(algword.coeffcients) == Set(coeffs)  #! test can be false if to equal monomoials are picked
                @test Set(algword.monomials) == Set(monom) broken = true
                @test ln == length(algword) #! test can be false if to equal monomoials are picked
            end
        end
    end
end
