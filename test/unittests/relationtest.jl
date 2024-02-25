import Nemo: CalciumField

k = CalciumField()
A, (a, b, c, d, e, f) = free_associative_algebra(k, ["a", "b", "c", "d", "e", "f"])

gen = VectorEnumeration.generators(A)

rel = [a^2*b + k(3)*e + k(1), (k(-1)/k(2))*d + k(-4) + k(7)*a + k(-9)*c, a^4*c^2*e - k(1), k(7)*a^3*b*c - (k(3)/k(4))*d*f^2 + a^2*d]

w = [VectorEnumeration.Weight(6,4), nothing, nothing, VectorEnumeration.Weight(10,5)]

fix = [true, true, false, true]

Relations = [VectorEnumeration.BinomialRelation(k(1), VectorEnumeration.GroupWord(a^2, gen), k(-3), VectorEnumeration.GroupWord(e, gen), VectorEnumeration.Weight(6,4)),
             VectorEnumeration.LinearRelation(VectorEnumeration.AlgebraWord((k(-1)/k(2))*d + k(-5) + k(7)*a + k(-9)*c, gen), VectorEnumeration.Weight(3)),
             VectorEnumeration.BinomialRelation(k(1), VectorEnumeration.GroupWord(a^4*c^2*e, gen), k(1), VectorEnumeration.GroupWord(A(1), gen), VectorEnumeration.Weight(3)),
             VectorEnumeration.AlgebraRelation(VectorEnumeration.AlgebraWord(k(7)*a^3*b*c - (k(3)/k(4))*d*f^2 + a^2*d, gen), VectorEnumeration.Weight(10,5))]



for i in 1:4
    tmp = VectorEnumeration.Relation(rel[i], gen; weight = w[i], fixing = fix[i])
    @test tmp == Relations[i] broken = true
end
