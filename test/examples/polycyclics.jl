@testset "Polycyclic groups" begin

    @testset "〈 x1, x2, x3 | x1^3 = x3, x2^2 = x3, x1°*x2*x1 = x2*x3, x1*x2*x1° = x2*x3 〉" begin
        #〈 x1, x2, x3 | x1^3*x3°, x2^2*x3°,x1°*x2*x1*x3°*x2°, x1*x2*x1°*x3°*x2° 〉

        field = QQ

        algebra, (x1, x2, x3, x1°, x2°, x3°) = free_associative_algebra(field, 6, Symbol("x"))

        inv = [x1 * x1° - 1, x1° * x1 - 1, x2 * x2° - 1, x2° * x2 - 1, x3 * x3° - 1, x3° * x3 - 1]

        rel = [x1^3 * x3° - 1, x2^2 * x3° - 1, x1° * x2 * x1 * x3° * x2° - 1, x1 * x2 * x1° * x3° * x2° - 1]

        relations = vcat(inv, rel)

        dimension = dimension_qa(algebra, relations)

        @test dimension == 6
    end

    @testset "ExamplesOfSomePcpGroups(1)" begin
        field = QQ

        algebra, (f1, f2, f3, f4, f1°, f2°, f3°, f4°) = free_associative_algebra(field, 8, Symbol("f"))

        inv = [f1 * f1° - 1, f1° * f1 - 1, f2 * f2° - 1, f2° * f2 - 1, f3 * f3° - 1, f3° * f3 - 1, f4 * f4° - 1, f4° * f4 - 1]

        rel = [f1° * f2 * f1 * f2° - 1, f1 * f2 * f1° * f2° - 1, f1° * f3 * f1 * f4° - 1, f1 * f3 * f1° * f4° * f3 - 1, f2° * f3 * f2 * f3 - 1, f2 * f3 * f2° * f3 - 1,
         f1° * f4 * f1 * f4° * f3° - 1, f1 * f4 * f1° * f3° - 1, f2° * f4 * f2 * f4 - 1, f2 * f4 * f2° * f4 - 1, f3° * f4 * f3 * f4° - 1, f3 * f4 * f3° * f4° - 1]

        relations = vcat(inv, rel)

        M = FreeModule(algebra, 1)

        subgens = [M([f1^10 - 1]), M([f2 - 1]), M([f3 - 1]), M([f4 - 1])] #Index = 10

        dimension = dimension_qm(algebra, M, relations, subgens)

        matrices = matrices_qm(SparseMatrixCSC, algebra, M, relations, subgens)
    
        # verify the relations are satisfied for low dimensions
        verify(relations, matrices)

        @testset "check dimension/order" begin
            @test dimension == 10
        end

        subgens = [M([f1^2 - 1]), M([f2 - 1]), M([f3 * f4^3 - 1]), M([f4^5 - 1])] #Index = 10

        dimension = dimension_qm(algebra, M, relations, subgens)

        matrices = matrices_qm(SparseMatrixCSC, algebra, M, relations, subgens)
    
        # verify the relations are satisfied for low dimensions
        verify(relations, matrices)

        @testset "check dimension/order" begin
            @test dimension == 10
        end
    end

    @testset "ExamplesOfSomePcpGroups(15)" begin
        field = QQ

        algebra, (f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13,
            f1°, f2°, f3°, f4°, f5°, f6°, f7°, f8°, f9°, f10°, f11°, f12°, f13°) = free_associative_algebra(field, 26, Symbol("f"))

        inv = [f1 * f1° - 1, f1° * f1 - 1, f2 * f2° - 1, f2° * f2 - 1, f3 * f3° - 1, f3° * f3 - 1, f4 * f4° - 1, f4° * f4 - 1,
        f5 * f5° - 1, f5° * f5 - 1, f6 * f6° - 1, f6° * f6 - 1, f7 * f7° - 1, f7° * f7 - 1, f8 * f8° - 1, f8° * f8 - 1, f9 * f9° - 1,
         f9° * f9 - 1, f10 * f10° - 1, f10° * f10 - 1, f11 * f11° - 1, f11° * f11 - 1, f12 * f12° - 1, f12° * f12 - 1, f13 * f13° - 1, f13° * f13 - 1]

        rel = [f1° * f2 * f1 * f3° * f2° - 1, f1 * f2 * f1° * f13°^5 * f11° * f10 * f9°^10 * f8^3 * f7° * f6° * f5 * f4° * f3 * f2° - 1, f1° * f3 * f1 * f4° * f3° - 1,
            f1 * f3 * f1° * f13 * f12°^2 * f11°^2 * f9°^3 * f6 * f5° * f4 * f3° - 1, f2° * f3 * f2 * f3° - 1, f2 * f3 * f2° * f3° - 1, f1° * f4 * f1 * f5° * f4° - 1, f1 * f4 * f1° * f6° * f5 * f4° - 1,
            f2° * f4 * f2 * f13°^4 * f11° * f10 * f9°^7 * f8^2 * f7° * f4° - 1, f2 * f4 * f2° * f13^6 * f12° * f11°^4 * f10°^2 * f9^7 * f8°^2 * f7 * f4° - 1, f3° * f4 * f3 * f13^6 * f12° * f11°^4 * f10°^2 * f9^7 * f8°^2 * f7 * f4° - 1,
            f3 * f4 * f3° * f13°^4 * f11° * f10 * f9°^7 * f8^2 * f7° * f4° - 1, f1° * f5 * f1 * f6° * f5° - 1, f1 * f5 * f1° * f6 * f5° - 1, f2° * f5 * f2 * f7° * f5° - 1, f2 * f5 * f2° * f13^6 * f12°^2 * f10°^3 * f7 * f5° - 1,
            f3° * f5 * f3 * f13 * f12°^3 * f11°^3 * f9°^5 * f8 * f5° - 1, f3 * f5 * f3° * f13^2 * f12°^2 * f11°^2 * f9^5 * f8° * f5° - 1, f4° * f5 * f4 * f13 * f12°^2 * f9^3 * f5° - 1, f4 * f5 * f4° * f13 * f12°^2 * f9°^3 * f5° - 1,
            f1° * f6 * f1 * f6° - 1, f1 * f6 * f1° * f6° - 1, f2° * f6 * f2 * f13^4 * f11°^4 * f10° * f9^7 * f8°^2 * f6° - 1, f2 * f6 * f2° * f12° * f11° * f10° * f9°^7 * f8^2 * f6° - 1, f3° * f6 * f3 * f13 * f12°^2 * f11°^3 * f9°^2 * f6° - 1,
            f3 * f6 * f3° * f13 * f12°^2 * f11°^2 * f9^2 * f6° - 1, f4° * f6 * f4 * f11°^2 * f6° - 1, f4 * f6 * f4° * f11°^3 * f6° - 1, f5° * f6 * f5 * f6° - 1, f5 * f6 * f5° * f6° - 1, f1° * f7 * f1 * f8° * f7° - 1,
            f1 * f7 * f1° * f11°^4 * f9° * f8 * f7° - 1, f2° * f7 * f2 * f13^6 * f12°^2 * f10°^3 * f7° - 1, f2 * f7 * f2° * f13°^4 * f12°^2 * f10^3 * f7° - 1, f3° * f7 * f3 * f13°^2 * f12° * f10 * f7° - 1, f3 * f7 * f3° * f13^4 * f12°^3 * f10° * f7° - 1,
            f4° * f7 * f4 * f13^2 * f12° * f7° - 1, f4 * f7 * f4° * f12°^3 * f7° - 1, f5° * f7 * f5 * f7° - 1, f5 * f7 * f5° * f7° - 1, f6° * f7 * f6 * f7° - 1, f6 * f7 * f6° * f7° - 1, f1° * f8 * f1 * f9° * f8° - 1,
            f1 * f8 * f1° * f11° * f9 * f8° - 1, f2° * f8 * f2 * f10° * f8° - 1, f2 * f8 * f2° * f10 * f8° - 1, f3° * f8 * f3 * f13 * f12°^3 * f8° - 1, f3 * f8 * f3° * f13 * f12° * f8° - 1, f4° * f8 * f4 * f8° - 1, f4 * f8 * f4° * f8° - 1, f5° * f8 * f5 * f8° - 1,
            f5 * f8 * f5° * f8° - 1, f6° * f8 * f6 * f8° - 1, f6 * f8 * f6° * f8° - 1, f7° * f8 * f7 * f8° - 1, f7 * f8 * f7° * f8° - 1, f1° * f9 * f1 * f11° * f9° - 1, f1 * f9 * f1° * f11°^4 * f9° - 1, f2° * f9 * f2 * f12° * f9° - 1,
            f2 * f9 * f2° * f13^2 * f12°^3 * f9° - 1, f3° * f9 * f3 * f9° - 1, f3 * f9 * f3° * f9° - 1, f4° * f9 * f4 * f9° - 1, f4 * f9 * f4° * f9° - 1, f5° * f9 * f5 * f9° - 1, f5 * f9 * f5° * f9° - 1, f6° * f9 * f6 * f9° - 1, f6 * f9 * f6° * f9° - 1,
            f7° * f9 * f7 * f9° - 1, f7 * f9 * f7° * f9° - 1, f8° * f9 * f8 * f9° - 1, f8 * f9 * f8° * f9° - 1, f1° * f10 * f1 * f13° * f10° - 1, f1 * f10 * f1° * f13 * f10° - 1, f2° * f10 * f2 * f10° - 1, f2 * f10 * f2° * f10° - 1,
            f3° * f10 * f3 * f10° - 1, f3 * f10 * f3° * f10° - 1, f4° * f10 * f4 * f10° - 1, f4 * f10 * f4° * f10° - 1, f5° * f10 * f5 * f10° - 1, f5 * f10 * f5° * f10° - 1, f6° * f10 * f6 * f10° - 1, f6 * f10 * f6° * f10° - 1, f7° * f10 * f7 * f10° - 1,
            f7 * f10 * f7° * f10° - 1, f8° * f10 * f8 * f10° - 1, f8 * f10 * f8° * f10° - 1, f9° * f10 * f9 * f10° - 1, f9 * f10 * f9° * f10° - 1, f11^5 - 1, f1° * f11 * f1 * f11° - 1, f1 * f11 * f1° * f11° - 1, f2° * f11 * f2 * f11° - 1, f2 * f11 * f2° * f11° - 1,
            f3° * f11 * f3 * f11° - 1, f3 * f11 * f3° * f11° - 1, f4° * f11 * f4 * f11° - 1, f4 * f11 * f4° * f11° - 1, f5° * f11 * f5 * f11° - 1, f5 * f11 * f5° * f11° - 1, f6° * f11 * f6 * f11° - 1, f6 * f11 * f6° * f11° - 1, f7° * f11 * f7 * f11° - 1,
            f7 * f11 * f7° * f11° - 1, f8° * f11 * f8 * f11° - 1, f8 * f11 * f8° * f11° - 1, f9° * f11 * f9 * f11° - 1, f9 * f11 * f9° * f11° - 1, f10° * f11 * f10 * f11° - 1, f10 * f11 * f10° * f11° - 1, f12^4 * f13°^2 - 1, f1° * f12 * f1 * f12° - 1,
            f1 * f12 * f1° * f12° - 1, f2° * f12 * f2 * f12° - 1, f2 * f12 * f2° * f12° - 1, f3° * f12 * f3 * f12° - 1, f3 * f12 * f3° * f12° - 1, f4° * f12 * f4 * f12° - 1, f4 * f12 * f4° * f12° - 1, f5° * f12 * f5 * f12° - 1, f5 * f12 * f5° * f12° - 1,
            f6° * f12 * f6 * f12° - 1, f6 * f12 * f6° * f12° - 1, f7° * f12 * f7 * f12° - 1, f7 * f12 * f7° * f12° - 1, f8° * f12 * f8 * f12° - 1, f8 * f12 * f8° * f12° - 1, f9° * f12 * f9 * f12° - 1, f9 * f12 * f9° * f12° - 1, f10° * f12 * f10 * f12° - 1,
            f10 * f12 * f10° * f12° - 1, f11° * f12 * f11 * f12° - 1, f1° * f13 * f1 * f13° - 1, f1 * f13 * f1° * f13° - 1, f2° * f13 * f2 * f13° - 1, f2 * f13 * f2° * f13° - 1, f3° * f13 * f3 * f13° - 1, f3 * f13 * f3° * f13° - 1,
            f4° * f13 * f4 * f13° - 1, f4 * f13 * f4° * f13° - 1, f5° * f13 * f5 * f13° - 1, f5 * f13 * f5° * f13° - 1, f6° * f13 * f6 * f13° - 1, f6 * f13 * f6° * f13° - 1, f7° * f13 * f7 * f13° - 1, f7 * f13 * f7° * f13° - 1, f8° * f13 * f8 * f13° - 1,
            f8 * f13 * f8° * f13° - 1, f9° * f13 * f9 * f13° - 1, f9 * f13 * f9° * f13° - 1, f10° * f13 * f10 * f13° - 1, f10 * f13 * f10° * f13° - 1, f11° * f13 * f11 * f13° - 1, f12° * f13 * f12 * f13° - 1]

        relations = vcat(inv, rel)

        F = FreeModule(algebra, 1)

        subgenlist = [[F([f1^20 - 1]), F([f2 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^4 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^8 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^12 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^16 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^2 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^6 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^10 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^14 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^18 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^5 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^9 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^13 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^17 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^3 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^7 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^11 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^15 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^19 - 1]), F([f2^20 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 - 1]), F([f2^10 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2^2 - 1]), F([f2^10 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2^4 - 1]), F([f2^10 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2^6 - 1]), F([f2^10 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2^8 - 1]), F([f2^10 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2 - 1]), F([f2^10 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2^3 - 1]), F([f2^10 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2^5 - 1]), F([f2^10 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2^7 - 1]), F([f2^10 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2^9 - 1]), F([f2^10 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^4 - 1]), F([f2^5 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^4 * f2 - 1]), F([f2^5 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^4 * f2^2 - 1]), F([f2^5 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^4 * f2^3 - 1]), F([f2^5 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^4 * f2^4 - 1]), F([f2^5 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^5 - 1]), F([f2^4 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^5 * f2^2 - 1]), F([f2^4 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^5 * f2 - 1]), F([f2^4 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^5 * f2^3 - 1]), F([f2^4 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^10 - 1]), F([f2^2 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^10 * f2 - 1]), F([f2^2 - 1]), F([f3 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^10 - 1]), F([f2 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^10 * f3 - 1]), F([f2 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 - 1]), F([f2^10 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 - 1]), F([f2^10 * f3 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^4 - 1]), F([f2^10 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])], [F([f1 * f2^4 - 1]), F([f2^10 * f3 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^8 - 1]), F([f2^10 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^8 - 1]), F([f2^10 * f3 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^2 - 1]), F([f2^10 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^2 - 1]), F([f2^10 * f3 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^6 - 1]), F([f2^10 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^6 - 1]), F([f2^10 * f3 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2 - 1]), F([f2^10 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2 - 1]), F([f2^10 * f3 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^5 - 1]), F([f2^10 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^5 - 1]), F([f2^10 * f3 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^9 - 1]), F([f2^10 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^9 - 1]), F([f2^10 * f3 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^3 - 1]), F([f2^10 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^3 - 1]), F([f2^10 * f3 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^7 - 1]), F([f2^10 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1 * f2^7 - 1]), F([f2^10 * f3 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 - 1]), F([f2^5 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f3 - 1]), F([f2^5 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2^2 - 1]), F([f2^5 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2^2 * f3 - 1]), F([f2^5 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2^4 - 1]), F([f2^5 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2^4 * f3 - 1]), F([f2^5 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2 - 1]), F([f2^5 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2 - 1]), F([f2^5 * f3 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2^3 - 1]), F([f2^5 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^2 * f2^3 - 1]), F([f2^5 * f3 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^5 - 1]), F([f2^2 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^5 - 1]), F([f2^2 * f3 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^5 * f2 - 1]), F([f2^2 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])],
            [F([f1^5 * f2 - 1]), F([f2^2 * f3 - 1]), F([f3^2 - 1]), F([f4 - 1]), F([f5 - 1]), F([f6 - 1]), F([f7 - 1]), F([f8 - 1]), F([f9 - 1]), F([f10 - 1]), F([f11 - 1]), F([f12 - 1]), F([f13 - 1])]]

        for subgens in subgenlist
            dimension = dimension_qm(algebra, F, relations, subgens)
            @test 20 == dimension
        end
    end
end




