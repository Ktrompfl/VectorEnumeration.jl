#Presentations and subgroup generators from GAP
@testset "Cosets" begin

    @testset "HS" begin
        #HS

        field = GF(17)

        algebra, (F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F1°, F3°, F4°, F5°, F6°, F8°, F10°) = free_associative_algebra(field, 17, Symbol("F"))

        inv = [F1*F1°-1, F1°*F1-1, F2^2-1, F3*F3°-1, F3°*F3-1, F4*F4°-1, F4°*F4-1, F5*F5°-1, F5°*F5-1, F6*F6°-1, F6°*F6-1, F7^2-1, F8*F8°-1, F8°*F8-1, F9^2-1, F10*F10°-1, F10°*F10-1]

        rel = [F1^3-1, F4°*F1*F3°*F2-1, (F2*F4)^2-1, F3^4-1, F4^2*F3*F2*F3°-1, F4*F3*F1*F4°*F3-1, (F5*F1°*F2)^2-1, F9*F2*F9*F3°*F4*F1-1, F6°*F3*F1*F5°*F2*F5-1, (F7*F3°)^3-1,
        F7*F2*F7*F4°*F3°*F1°-1, (F5°*F3)^3-1, F7*F1°*F3°*F7*F3*F1-1, F5*F1°*F5^2*F2*F5°*F1°-1, F7*F1*F4°*F7*F1°*F4*F1-1, F9*F3°*F9*F4°*F2*F6*F1°*F4-1,
        F5*F7*F2*F7*F5*F6°*F1°*F3°-1, F7*F3°*F7*F3*F7*F1*F3*F1-1, F5*F1°*F3*F6*F5°^2*F6°*F4-1, F8°*F5*F3*F5*F1*F7*F5°*F3°-1, F9*F5°*F9*F1*F2*F5°*F3^2*F1°-1,
        F5*F1*F2*F5*F1*F5°*F3°*F1°*F4°-1, F9*F1°*F9*F5*F1°*F5*F3°*F1*F4*F2-1, F10*F3°*F10°*F1°*F5*F3*F5*F3°*F1*F2-1, F7*F5°*F4*F7*F5*F7*F3*F4*F5°*F3*F5°-1,
        F3*F5°*F1*F5°*F7*F5*F1°*F5*F7*F4°*F1-1, F10*F7*F10°*F3*F5*F7*F4°*F1*F5*F1°*F4°*F1°-1, F10*F8°*F5*F10*F1*F6*F5°*F7*F5°*F4*F3°*F1-1,
        F10*F1°*F10°*F5^2*F7*F5*F4°*F5*F1*F4*F1-1, F10°*F2*F9*F10*F4*F1°*F5*F4*F5°*F7*F5*F4*F5°-1, F10*F2*F10°*(F5°*F7)^2*F2*F5°*F3°*F5*F4*F3°*F1°-1,
        F10*F5°*F8°*F7*F9*F4*F3°*F5*F3*F2*F7*F5°*F4*F5°-1, F10*F6°*F1°*F7*F9*(F4*F5°)^2*F1°*F7*F5°*F4*F5°-1, F10°*F4*F2*F10*F3°*F10°*F4°*F2*F10*F1°*F3°*F1*F3^2*F1-1,
        F10°*F4°*F2*F10*F3*F5°*F4°*F5°*F7*F5°*F3°*F5*F3^2*F1-1, F10*F6*F9*F10°*F2*F10*F3°*F10°*F3°*F4*F7*F5°*F4°*F5°*F7-1,
        F10*F6*F1*F7*F9*F3*F5°*F4°*F5°*F7*(F5*F4°)^2*F3-1, F10°*F6°*F5*F8°*F2*F9*F3^2*F5*F4*F7*F5°^2*F1°*F5°-1] #Group order = 44352000

        relations = vcat(inv, rel)

        F = FreeModule(algebra, 1)

        @testset "Quotient 1" begin
            subgens = [F([F1°*F4°*F1*F5°*F6*F5°*F8^2*F2*F10*F8*F5°-1]), F([F1°*F3*F4°*F3°*F5*F1*F7*F5*F3*F10*F8*F6°-1])] #Index = 100
            dimension = dimension_qm(algebra, F, relations, subgens)

            @testset "check dimension/order" begin
                @test dimension == 100
            end

        end
        @testset "Quotient 2" begin

            subgens = [F([F1°*F3°^2*F1*F5^2*F1°*F8*F2-1]), F([F3*F4°*F1°*F10*F3°*F10°*F8°*F9*F8*F1°-1])] #Index = 100
            dimension = dimension_qm(algebra, F, relations, subgens)

            @testset "check dimension/order" begin
                @test dimension == 100
            end

        end
        @testset "Quotient 3" begin

            subgens = [F([F1*F4°*F5*F3*F5°*F8*F1*F2°*F10^2*F7*F1-1]), F([F4°*F3°*F6°*F5°*F8*F4*F1°-1])] #Index = 176
            dimension = dimension_qm(algebra, F, relations, subgens)

            @testset "check dimension/order" begin
                @test dimension == 176
            end

        end
        @testset "Quotient 4" begin

            subgens = [F([F1*F2°*F5°*F6*F5°*F1°*F7*F8°*F9*F8*F3°-1]), F([F1*F3*F4°*F5°*F6*F1°*F8°*F5°*F4*F10*F5*F3-1])] #Index = 1100
            dimension = dimension_qm(algebra, F, relations, subgens)

            @testset "check dimension/order" begin
                @test dimension == 1100
            end

        end
    end

    @testset "A₁₀" begin


        #A₁₀


        field = QQ

        algebra, (A_101, A_102, A_101°, A_102°) = free_associative_algebra(field, 4, Symbol("A10"))

        inv = [A_101*A_101°-1, A_101°*A_101-1, A_102*A_102°-1, A_102°*A_102-1]

        rel = [A_101^9*A_102°^9-1, A_101^9*(A_102°*A_101°)^5-1, (A_101°*A_102°*A_101*A_102)^2-1, (A_101°^2*A_102°*A_101*A_102^2)^2-1, (A_101°^3*A_102°*A_101*A_102^3)^2-1,
            (A_101°^4*A_102°*A_101*A_102^4)^2-1] # Order = 1814400

        relations = vcat(inv, rel)

        F = FreeModule(algebra, 1)

        @testset "Quotient 1" begin

            subgens = [F([A_102*A_101*A_102*A_101°^4*A_102*A_101-1]), F([A_101^2*A_102^2*A_101°*A_102°*A_101°^2-1]), F([A_101*A_102°-1]), F([A_101°*A_102-1]),
                F([A_102*(A_101°*A_102°)^2*A_101*A_102^2*A_101*A_102°-1])] #Index = 210

            dimension = dimension_qm(algebra, F, relations, subgens)

            @testset "check dimension/order" begin
                @test dimension == 210 
            end
        end

        @testset "Quotient 2" begin

            subgens = [F([A_102*A_101*A_102°*A_101°*(A_101°*A_102°)^2*(A_101*A_102)^2*A_101-1]), F([A_102°*A_101^3*A_102°^2-1]), F([A_101°^4*A_102°^4*A_101°-1]),
                F([A_101*A_102°^2*A_101°*A_102°*A_101*A_102*A_101-1])] #Index = 120

            dimension = dimension_qm(algebra, F, relations, subgens)

            @testset "check dimension/order" begin
                @test dimension == 120 
            end
        end
    end

    #TODO: adjust to new functions (base_qm, ...)
    #TODO: currently not terminating
    #=
    @testset "GL(4, 3)" begin

        #GL(4, 3)


        field = QQ

        algebra, (F1 , F2 , F3 , F4 , F5 , F6 , F7 , F8 , F9 , F10 , F11 , F12 , F13 , F14 , F15 , F16,
                F1°, F2°, F3°, F4°, F5°, F6°, F7°, F8°, F9°, F10°, F11°, F12°, F13°, F14°, F15°, F16°  ) = free_associative_algebra(field, 32, Symbol("F"))

        T = elem_type(algebra)

        inv = Union{T, Tuple{T,T}}[(F1,F1°), (F2,F2°), (F3,F3°), (F4,F4°), (F5,F5°), (F6,F6°), (F7,F7°), (F8,F8°), (F9,F9°), (F10,F10°), (F11, F11°), (F12, F12°), (F13, F13°),
        (F14, F14°), (F15, F15°), (F16, F16°)]

        rel = Union{T, Tuple{T,Weight}}[ F1^2, F1°*F2*F1*F3*F16°, F1°*F3*F1*F2*F16°, F1°*F4*F1*F4°, F1°*F5*F1*F6°*F5*F16°, F1°*F6*F1*F6, F1°*F7*F1*F4°*F3°*F5, 
        F1°*F8*F1*F11*F8*F11°*F2*F9°*F7°*F6*F16°, F1°*F9*F1*F5°*F2*F4*F9*F12°*F6°, F1°*F10*F1*F15*F14°*F3*F6°*F5*F16°, F1°*F11*F1*F12*F2*F4°*F9°*F2*F6, 
        F1°*F12*F1*F5°*F6°*F12°*F6°, F1°*F13*F1*F13°*F12*F10°*F9*F5*F16°, F1°*F14*F1*F14°, F1°*F15*F1*F5°*F10*F2°*F14*F16°, F2^2, F3^2, F5^2, F7^2, F8^2, 
        F12^2*F9*F16°, F9^3*F16°, (F14°*F2)^2, F8*F9*F8*F9°, (F9°*F5)^2, (F7*F4)^2, (F13°*F9)^2, (F5*F3)^2, (F7*F5)^2, (F2*F3)^2, F9*F7*F9°*F7, F14*F12*F9°*F15°*F2, 
        F6°*F7*F5*F4°*F2, F4°*F7*F4*F5*F3, F10°*F6°*F9°*F2*F4*F16°, F8*(F9*F3)^2*F7, F4*F5*F4*F3*F4°*F3, (F8*F2*F4°)^2, F8*F3*F8*F2*F7*F5*F2, F5*(F2*F4°)^2*F2*F3*F16°, 
        F8*F9°*F4°*F9*F3*F4*F3, F9°*F6°*F8*F2*F9°*F2*F6°, F9*F4*F10°*F4°*F9°*F2*F7*F3*F16°, F11*F2*F4°*F11°*F4°*F6°*F9*F4*F16°, F2*F15°*F4*F6*F11*F6*F9*F4*F14°
        , F11*F3*F11°*F9*F2*F3*F9*F2*F5, F14*F15*F4*F5*F11*F6*F9*F4*F3, F11*F9°*F2*F11*F5*F4*F6*F9°*F2*F16°, F13*F2*F10*F13*F9*F11*F3*F6*F9, F12°*F4*F5*F11*F2*F9*F2*F4°*F5*F16°,
        F13°*F4°*F9*F8*F13°*F8*F11*F6°*F4*F3*F16°, (F13°*F15)^2*F2*F14*F4°*F6°*F2*F5*F16°, F11°*F4*F7*F11*F2*F8*F3*F7*F4°*F5*F16°, 
        F13°*F3*F4*F14°*F9*F11*F6°*F9°*F6°*F4°, F13*F6°*F13°*F5*F4*F9*F6*F11°*F9°*F2*F16°, F13*F4*F5*F13°*F3*F4*F9°*F6°*F11°*F6*F16°, 
        F13°*F6°*F14*F13*F4°*F6*F11*F6*F9*F4°*F16°, F13°*F6°*F14°*F13*F4°*F11*F6*F9*F6*F7*F16°, F11*F2*F4*F11°*F2*F6*F9*F4°*F3*F7, 
        F13*F14*F15*F13*F11*F6*F9*F6°*F3*F7*F16°, F11^2*F9*F11*F3*F4°*F9*F6*F2*F7*F5, F13*F4^2*F13°*F5*F3*F11*F6*F9*F7*F5, F14°*F11°*F13°*F5*F4*F2*F9*F2*F4°*F11*F8*F9*F16°, 
        F13*F11°*F2*F13*F5*F2*F9*F7*F2*F4°*F11°*F9, F13*F3*F15°*F13°*F9°*F11*F9*F6*F9*F2*F3*F6°, F13*F10°*F4°*F13*F11*F2*F11°*F2*F6*F9°*F6°*F3, 
        F13*F10°*F12°*F13°*F4*F3*F11*F3*F2*F9*F4*F2*F7*F5*F16°, F1°*F16*F1*F16°, F2°*F16*F2*F16°, F3°*F16*F3*F16°, F4°*F16*F4*F16°, F5°*F16*F5*F16°, 
        F6°*F16*F6*F16°, F7°*F16*F7*F16°, F8°*F16*F8*F16°, F9°*F16*F9*F16°, F10°*F16*F10*F16°, F11°*F16*F11*F16°, F12°*F16*F12*F16°, F13°*F16*F13*F16°, 
        F14°*F16*F14*F16°, F15°*F16*F15*F16°, F16^2 ]

        @testset "Quotient 1" begin

            subgens = NTuple{1, T}[ (F5°*F2*F4°*F10°*F9*F8°*F11°*F13-1,), (F6°*F14*F9*F6°*F12*F10*F4*F7-1,), (F16*F8*F11°*F4°*F9*F7°*F2°*F1-1,), (F16*F6*F3°*F10°*F9*F11°*F4°*F13*F8°-1,) ] #Index = 130
            (generators, matrices), (dimension, ntotal), (img, preimg) = vector_enumeration(algebra, inv, rel, 1, subgens)

            verify(field, inv, rel, dimension, matrices)

            @testset "check dimension/order" begin            
                @test dimension == 100
            end

        end

    end

    =#
end
