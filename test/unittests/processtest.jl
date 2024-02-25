k = GF(23)
A, (a, b, c, d) = free_associative_algebra(k, ["a", "b", "c", "d"])
gen = VectorEnumeration.generators(A)



rel = [ VectorEnumeration.Relation(a*b^2, gen),
        VectorEnumeration.Relation(c^2, gen), 
        VectorEnumeration.Relation(b + k(2) * c, gen), 
        VectorEnumeration.Relation(b^3 - k(4)a*c, gen), 
        VectorEnumeration.Relation(d^2 - (a*c)^2 + 1, gen)]

sub = NTuple{2, Union{VectorEnumeration.AlgebraWord, VectorEnumeration.GroupWord, VectorEnumeration.Generator}}[ (VectorEnumeration.AlgebraWord(a^2 + k(2), gen), VectorEnumeration.GroupWord(one(A), gen))]

a, b, c, d = gen

c.inverse = c

basis = VectorEnumeration.Basis(k)

start_basis = [ VectorEnumeration.create!(basis ,VectorEnumeration.BasisElementData(1, VectorEnumeration.GroupWord(one(A), gen), 4, 1)),
                VectorEnumeration.create!(basis ,VectorEnumeration.BasisElementData(2, VectorEnumeration.GroupWord(one(A), gen), 4, 1))]

b₁, b₂ = start_basis
z = zero(b₁)



@testset "process_subgen!" begin
    VectorEnumeration.process_subgen!(sub[1], start_basis, 1, gen)
    @test iszero(b₁[a][a] + k(2)*b₁ + b₂)
    @test basis.created == 4
    @test basis.deleted == 1
    @test VectorEnumeration.isdeleted(basis.last)
end

@testset "define, replace, action, process_equation" begin
    #processing relation b₁.a*b^2 = b₁ 
    v = VectorEnumeration.action!(b₁, a, 2, gen) #returning b₁[a] as already defined
    b₃ = b₁[a]
    @test basis.created == 4
    @test iszero(v - b₁[a])
    VectorEnumeration.define_action!(v.elem[1], b, 2, gen)
    @test basis.created == 5
    @test iszero(v.elem[1][b] - basis.last)
    v = VectorEnumeration.action!(v, b, 2, gen)
    b₅ = v.elem[1]
    @test basis.created == 5
    v = VectorEnumeration.action!(v, b, 2, gen)
    @test basis.created == 6
    b₆ = v.elem[1] 
    VectorEnumeration.replace!(b₆, b₁, 2, gen) # v = b₆ =ꜝ b₁
    @test VectorEnumeration.isdeleted(b₆)
    @test basis.created == 6
    @test basis.deleted == 2
    @test VectorEnumeration.replacement(b₆) === b₁
    @test VectorEnumeration.iszero(VectorEnumeration.action!(b₁, a*b*b, 2, gen) - b₁)
    #process b₁.c^2 = b₁
    VectorEnumeration.define_action!(b₁, c, 2, gen)
    @test basis.created == 7
    b₇ = b₁[c]
    @test b₇ == basis.last == basis.last_undeleted
    @test !VectorEnumeration.isnothing(b₇[c])
    v = VectorEnumeration.action!(b₁, c*c, 2, gen)
    @test basis.created == 7
    @test basis.deleted == 2
    VectorEnumeration.process_equation!(v, b₁, 2, gen) 
    @test VectorEnumeration.iszero(v - b₁)
    @test basis.deleted == 2
    rel[1].processed_last = b₁
    rel[2].processed_last = b₁
end
@testset "action, process_equation, process_generator_equation" begin
    #process b₁.(b + 2c) = b₁
    v = VectorEnumeration.action!(b₁, rel[3].word, 3, gen)
    @test basis.created == 8
    VectorEnumeration.process_equation!(v, z, 3, gen)
    @test basis.created == 8
    @test basis.deleted == 3
    @test VectorEnumeration.iszero(b₁[b] + k(2)*b₁[c] - b₁)
    rel[3].processed_last = b₁
    #process b₁.(b^3 - 4ac) = b₁
    v = VectorEnumeration.action!(b₁, rel[4].word, 3, gen)
    @test basis.created == 11
    b₁₁ = v.elem[1]
    VectorEnumeration.define_action!(b₁₁, d, 3, gen)
    #force process generator equation by defining b₁₁[d]
    @test basis.created == 12
    VectorEnumeration.process_equation!(v, b₁, 3, gen)
    @test basis.created == 15
    @test basis.deleted == 4
    @test VectorEnumeration.iszero(v - b₁)
    @test VectorEnumeration.iszero(VectorEnumeration.action!(v - b₁, d, 3, gen))
    @test basis.created == 15
    rel[4].processed_last = b₁
end

@testset "process_relation" begin
    VectorEnumeration.process_relation!(b₁, rel[5], 3, gen)
    @test VectorEnumeration.iszero(VectorEnumeration.action!(rel[5].coeff₁*b₁, rel[5].word₁, 3, gen)-VectorEnumeration.action!(rel[5].coeff₂*b₁, rel[5].word₂, 3, gen))
    rel[5].processed_last = b₁
    b = VectorEnumeration.next(b₁)
    VectorEnumeration.process_relation!(b, rel[1], 4, gen)
    @test VectorEnumeration.iszero(VectorEnumeration.action!(b, rel[1].word₁, 4, gen)-b)
    VectorEnumeration.process_relation!(b, rel[2], 4, gen)
    @test VectorEnumeration.iszero(VectorEnumeration.action!(b, rel[2].word₁, 4, gen)-b)
    VectorEnumeration.process_relation!(b, rel[3], 4, gen)
    @test VectorEnumeration.iszero(VectorEnumeration.action!(b, rel[3].word, 4, gen))
    println(VectorEnumeration.isdeleted(b))
    VectorEnumeration.process_relation!(b, rel[4], 4, gen)
    @test VectorEnumeration.iszero(VectorEnumeration.action!(b, rel[4].word, 4, gen)-b)
end

