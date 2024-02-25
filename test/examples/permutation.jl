#= permutation representations of various small groups =#

### symmetric groups
N = 7
@testset "symmetric groups Sₙ for 3 ≤ n ≤ $N" begin
    for n = 3:N
        @testset "symmetric group Sₙ⁽¹⁾ for n = $n" begin
            # presentation taken from Coxeter Graph Aₙ in https://agag-lassueur.math.rptu.de/~lassueur/en/teaching/ProSemSS22/Coxeter_WS1920.pdf
            A, X = free_associative_algebra(QQ, n - 1, "x")
            inv = [X[i]^2- 1 for i = 1:n-1]
            rel = []
            for i in 1:(n-2)
                push!(rel, (X[i] * X[i+1])^3 -1)
                for j in (i+2):(n-1)
                    push!(rel, (X[i] * X[j])^2 -1)
                end
            end

            relations = Vector{elem_type(A)}(vcat(inv, rel))

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(SparseMatrixCSC, A, relations)

            verify(relations, matrices)

            # test dimension, i.e. order of group
            @test dimension == factorial(n)  # |Sₙ| = n!
        end

        @testset "symmetric group Sₙ⁽²⁾ for n = $n" begin
            # presentation from https://math.stackexchange.com/questions/3972026/what-are-the-relations-in-this-presentation-of-s-n
            #⟨t, c | t² = cⁿ = (tc)ⁿ⁻¹ = [t,c]³ = 1; [t,cᵏ]² = 1 for 2 ≤ k ≤ n÷2⟩
            A, X = free_associative_algebra(QQ, 3, "x")
            t = X[1]
            c = X[2]
            c⁻¹ = X[3]
            inv = [t^2 - 1, c*c⁻¹ - 1, c⁻¹*c - 1]
            rel = [c^n - 1, (t * c)^(n - 1) - 1, (t * c⁻¹ * t * c)^3 - 1]
            for k in 2:n÷2
                push!(rel, (t * c⁻¹^k * t * c^k)^2 - 1)
            end

            relations = Vector{elem_type(A)}(vcat(inv, rel))

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(SparseMatrixCSC, A, relations)

            verify(relations, matrices)

            # test dimension, i.e. order of group
            @test dimension == factorial(n)  # |Sₙ| = n!
        end
    end
end


### alternating groups
N = 8
@testset "alternating groups Aₙ for 3 ≤ n ≤ $N" begin
    for n = 3:N
        @testset "alternating group Aₙ for n = $n" begin
            # presentation from https://math.stackexchange.com/questions/152158/presentation-of-a-n-from-jacobsons-basic-algebra-i
            A, X = free_associative_algebra(QQ, (n - 1), "x")
            inv = [X[1]*X[n-1] - 1, X[n-1]*X[1] - 1]
            rel = [X[1]^3 - 1]
            for i in 1:(n-3)
                push!(inv, X[i+1]^2 - 1)
                push!(rel, (X[i] * X[i+1])^3 - 1)
                for j in (i+2):(n-2)
                    push!(rel, (X[i] * X[j])^2 - 1)
                end
            end

            relations = Vector{elem_type(A)}(vcat(inv, rel))

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(SparseMatrixCSC, A, relations)

            # verify the relations are satisfied for low dimensions
            verify(relations, matrices)

            # test dimension, i.e. order of group
            @test dimension == factorial(n) ÷ 2  # |Aₙ| = n!/2
        end
    end
end


# dihedral groups
N = 10
@testset "dihedral groups D₂ₙ for 2 ≤ n ≤ $N" begin
    for n =2:N 
        @testset "dihedral group D₂ₙ⁽¹⁾ for n = $n" begin
            # ⟨ a,b | a^n, b^2, (a*b)^2 ⟩ = D₂ₙ
            A, X = free_associative_algebra(QQ, 3, "x")
            a = X[1]
            b = X[3]
            inv = [X[1]*X[2]-1, X[2]*X[1]-1, b^2-1]
            rel = [a^n - 1, (a * b)^2 - 1]

            relations = Vector{elem_type(A)}(vcat(inv, rel))

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(SparseMatrixCSC, A, relations)

            # verify the relations are satisfied for low dimensions
            verify(relations, matrices)

            # test dimension, i.e. order of group
            @test dimension == 2 * n  # |D₂ₙ| = 2n
        end

        @testset "dihedral group D₂ₙ⁽²⁾ for n = $n" begin
            # ⟨ a,b | a^2, b^2, (a*b)^n ⟩ = D₂ₙ
            A, X = free_associative_algebra(QQ, 2, "x")
            a = X[1]
            b = X[2]
            inv = [a^2 - 1, b^2 - 1]
            rel = [(a * b)^n - 1]

            relations = Vector{elem_type(A)}(vcat(inv, rel))

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(SparseMatrixCSC, A, relations)

            # verify the relations are satisfied for low dimensions
            verify(relations, matrices)

            # test dimension, i.e. order of group
            @test dimension == 2 * n  # |D₂ₙ| = 2n
        end
    end
end
 

### some of the presentations used to compare different implementations of coset enumeration in 
# Cannon, John J., Lucien A. Dimino, George Havas, and Jane M. Watson. 
# “Implementation and Analysis of the Todd-Coxeter Algorithm.” 
# Mathematics of Computation 27, no. 123 (1973): 463–90. 
# https://doi.org/10.2307/2005654.
𝔽₂ = QQ # TODO change to F₂ #?

@testset "E₁" begin
    A, (s, t, r, s′, t′, r′) = free_associative_algebra(𝔽₂, ["s", "t", "r", "s'", "t'", "r'"])
    inv = [s*s′ - 1, t*t′ - 1, r*r′ - 1, s′*s - 1, t′*t - 1, r′*r - 1]
    rel = [t′ * r * t * r′^2 - 1, r′ * s * r * s′^2 - 1, s′ * t * s * t′^2 - 1]

    relations = Vector{elem_type(A)}(vcat(inv, rel))

    dimension = dimension_qa(A, relations)

    matrices = matrices_qa(SparseMatrixCSC, A, relations)

    # verify the relations are satisfied for low dimensions
    verify(relations, matrices)

    # test dimension, i.e. order of group
    @test dimension == 1
end

@testset "Cox" begin
    A, (a, b, a′, b′) = free_associative_algebra(𝔽₂, ["a", "b", "a'", "b'"])
    inv = [a*a′ - 1, b*b′ - 1, a′*a - 1, b′*b - 1]
    rel = [a^6 - 1, b^6 - 1, (a * b)^2 - 1, (a^2 * b^2)^2 - 1, (a^3 * b^3)^5 - 1]

    relations = Vector{elem_type(A)}(vcat(inv, rel))

    dimension = dimension_qa(A, relations)

    matrices = matrices_qa(SparseMatrixCSC, A, relations)

    # verify the relations are satisfied for low dimensions
    #verify(relations, matrices)

    # test dimension, i.e. order of group
    @test dimension == 3000
end

@testset "B₂,₄" begin
    A, (a, b, a′, b′) = free_associative_algebra(𝔽₂, ["a", "b", "a'", "b'"])
    inv = [a*a′ - 1, b*b′ - 1, a′*a - 1, b′*b - 1]
    rel = [a^4 - 1, b^4 - 1, (a * b)^4 - 1, (a′ * b)^4 - 1, (a^2 * b)^4 - 1, (a * b^2)^4 - 1, (a^2 * b^2)^4 - 1, (a′ * b * a * b)^4 - 1, (a * b′ * a * b)^4 - 1]

    relations = Vector{elem_type(A)}(vcat(inv, rel))

    dimension = dimension_qa(A, relations)

    matrices = matrices_qa(SparseMatrixCSC, A, relations)

    # verify the relations are satisfied for low dimensions
    #verify(relations, matrices)

    # test dimension, i.e. order of group
    @test dimension == 4096
end

@testset "S₇" begin
    A, (a, b, a′, b′) = free_associative_algebra(𝔽₂, ["a", "b", "a'", "b'"])
    inv = [a*a′ - 1, b*b′ - 1, a′*a - 1, b′*b - 1]
    rel = [a^7 - 1, b^2 - 1, (a * b)^6 - 1, (a * b * a′ * b′)^3 - 1, (a * a * b * a′ * a′ * b′)^2 - 1, (a * a * a * b * a′ * a′ * a′ * b′)^2 - 1]

    relations = Vector{elem_type(A)}(vcat(inv, rel))

    dimension = dimension_qa(A, relations)

    matrices = matrices_qa(SparseMatrixCSC, A, relations)

    # verify the relations are satisfied for low dimensions
    verify(relations, matrices)

    # test dimension, i.e. order of group
    @test dimension == 5040
end


@testset "PSL₂(11)" begin
    A, (a, b, a′) = free_associative_algebra(𝔽₂, ["a", "b", "a'"])
    inv = [a*a′ - 1, a′*a - 1, b^2 - 1]
    rel = [a^11-1, (a * b)^3-1, (a^4 * b * a′^5 * b)^2-1]

    relations = Vector{elem_type(A)}(vcat(inv, rel))

    dimension = dimension_qa(A, relations)

    matrices = matrices_qa(SparseMatrixCSC, A, relations)

    # verify the relations are satisfied for low dimensions
    #verify(relations, matrices)

    # test dimension, i.e. order of group
    @test dimension == 660
end

@testset "PSL₂(17)" begin
    A, (a, b, a′) = free_associative_algebra(𝔽₂, ["a", "b", "a'"])
    inv = [a*a′ - 1, a′*a - 1, b^2 - 1]
    rel = [a^9-1, (a * b)^4-1, (a^2 * b)^3-1]

    relations = Vector{elem_type(A)}(vcat(inv, rel))

    dimension = dimension_qa(A, relations)

    matrices = matrices_qa(SparseMatrixCSC, A, relations)

    # verify the relations are satisfied for low dimensions
    #verify(relations, matrices)

    # test dimension, i.e. order of group
    @test dimension == 2448
end

@testset "PSL₃(4)" begin
    A, (a, b, a′, b′) = free_associative_algebra(𝔽₂, ["a", "b", "a'", "b'"])
    inv = [a*a′ - 1, b*b′ - 1, a′*a - 1, b′*b - 1]
    rel = [a^5-1, b^3 - 1, (a * b)^4 - 1, (a′ * b * a′ * b′ * a * b)^3 - 1, (b * a^2 * b * a′^2 * b′ * a′^2 * b′ * a^2 * b * a * b′ * a′) - 1]

    relations = Vector{elem_type(A)}(vcat(inv, rel))

    dimension = dimension_qa(A, relations)

    matrices = matrices_qa(SparseMatrixCSC, A, relations)

    # verify the relations are satisfied for low dimensions
    #verify(relations, matrices)

    # test dimension, i.e. order of group
    @test dimension == 20160
end

@testset "M₁₁⁽¹⁾" begin
    A, (a, b, c, a′, b′, c′) = free_associative_algebra(𝔽₂, ["a", "b", "c", "a'", "b'", "c'"])
    inv = [a*a′ - 1, b*b′ - 1, c*c′ - 1, a′*a - 1, b′*b - 1, c′*c - 1]
    rel = [a^11 - 1, b^5 - 1, c^4 - 1, (a^4 * c^2)^3 - 1, (b * c^2)^2 - 1, (a * b * c)^3 - 1, b′ * a * b * a′^4 - 1, c′ * b * c * b′^2 - 1]

    relations = Vector{elem_type(A)}(vcat(inv, rel))

    dimension = dimension_qa(A, relations)

    matrices = matrices_qa(SparseMatrixCSC, A, relations)

    # verify the relations are satisfied for low dimensions
    #verify(relations, matrices)

    # test dimension, i.e. order of group
    @test dimension == 7920
end

@testset "M₁₁⁽²⁾" begin
    A, (a, b, c, a′, b′, c′) = free_associative_algebra(𝔽₂, ["a", "b", "c", "a'", "b'", "c'"])
    inv = [a*a′ - 1, b*b′ - 1, c*c′ - 1, a′*a - 1, b′*b - 1, c′*c - 1]
    rel = [a^11 - 1, b^5 - 1, c^4 - 1, (a * c)^3 - 1, c′ * b * c * b′^2 - 1, b′ * a * b * a′^4 - 1]

    relations = Vector{elem_type(A)}(vcat(inv, rel))

    dimension = dimension_qa(A, relations)

    matrices = matrices_qa(SparseMatrixCSC, A, relations)

    # verify the relations are satisfied for low dimensions
    #verify(𝔽₂, inv, rel, dimension, matrices)

    # test dimension, i.e. order of group
    @test dimension == 7920
end

# this test takes quite long to complete
# @testset "Neu" begin
#     A, (a, b, c, a′, b′, c′) = free_associative_algebra(𝔽₂, ["a", "b", "c", "a'", "b'", "c'"])
#     inv = [a*a′ - 1, b*b′ - 1, c*c′ - 1, a′*a - 1, b′*b - 1, c′*c - 1]
#     rel = [ a^3- 1, b^3 - 1, c^3 - 1, (a*b)^5 - 1, (a′*b)^5 - 1, (a*c)^4 - 1, (a*c′)^4 - 1, a*b′*a*b*c′*a*c*a*c′ - 1, (b*c)^3 - 1, (b′*c)^4 - 1]

#     relations = Vector{elem_type(A)}(vcat(inv, rel))

#     dimension = dimension_qa(A, relations)

#     matrices = matrices_qa(SparseMatrixCSC, A, relations)

#     # verify the relations are satisfied for low dimensions
#     verify(𝔽₂, inv, rel, dimension, matrices)

#     # test dimension, i.e. order of group
#     @test dimension == 40320
# end

@testset "Weyl B₆" begin
    A, (a, b, c, d, e, f) = free_associative_algebra(𝔽₂, ["a", "b", "c", "d", "e", "f"])
    inv = [a^2 - 1, b^2 - 1, c^2 - 1, d^2 - 1, e^2 - 1, f^2 - 1]
    rel = [(a * b)^3 - 1, (a * c)^2 - 1, (a * d)^2 - 1, (a * e)^2 - 1, (a * f)^2 - 1, (b * c)^3 - 1, (b * d)^2 - 1, (b * e)^2 - 1, (b * f)^2 - 1, (c * d)^3 - 1, (c * e)^2 - 1, (c * f)^2 - 1, (d * e)^3 - 1, (d * f)^2 - 1, (e * f)^4 - 1]

    relations = Vector{elem_type(A)}(vcat(inv, rel))

    dimension = dimension_qa(A, relations)

    matrices = matrices_qa(SparseMatrixCSC, A, relations)

    # verify the relations are satisfied for low dimensions
    #verify(𝔽₂, inv, rel, dimension, matrices)

    # test dimension, i.e. order of group
    @test dimension == 46080
end



### three polynomial presentations of the form ⟨s, t, u | stu - 1, p₁(s), p₂(t), p₃(u)⟩ for polynomials p₁, p₂ and p₃ from  
# S.A. Linton, Constructing matrix representations of finitely presented groups,
# Journal of Symbolic Computation, Volume 12, Issues 4–5, 1991,Pages 427-438, ISSN 0747-7171,
# https://doi.org/10.1016/S0747-7171(08)80095-8. (https://www.sciencedirect.com/science/article/pii/S0747717108800958)
@testset "three polynomial presentations" begin
    @testset "examples over GF(2)" begin
        𝔽₂ = GF(2)
        A, (s, t, u, s⁻¹, t⁻¹, u⁻¹) = free_associative_algebra(𝔽₂, ["s", "t", "u", "s⁻¹", "t⁻¹", "u⁻¹"])
        inv = [s * s⁻¹ - 1, s⁻¹ * s - 1, t * t⁻¹ - 1, t⁻¹ * t - 1, u * u⁻¹ - 1, u⁻¹ * u - 1]

        @testset "p₁ = x² + 1, p₂ = x² + x + 1, p₃ = x⁷ + 1" begin
            rel = [s * t * u + 1, s^2 + 1, t^2 + t + 1, u^7 + 1]

            relations = vcat(inv, rel)

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(Matrix, A, relations)

            # verify the relations are satisfied for low dimensions
            verify(relations, matrices)

            # test dimension
            @test dimension == 12
        end

        @testset "p₁ = x² + 1, p₂ = x² + x + 1, p₃ = x⁹ + 1" begin
            rel = [s * t * u + 1, s^2 + 1, t^2 + t + 1, u^9 + 1]

            relations = vcat(inv, rel)

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(Matrix, A, relations)

            # verify the relations are satisfied for low dimensions
            verify(relations, matrices)

            # test dimension
            @test dimension == 16
        end

        @testset "p₁ = x² + 1, p₂ = x² + x + 1, p₃ = x¹¹ + 1" begin
            rel = [s * t * u + 1, s^2 + 1, t^2 + t + 1, u^11 + 1]

            relations = vcat(inv, rel)

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(Matrix, A, relations)

            # verify the relations are satisfied for low dimensions
            verify(relations, matrices)

            # test dimension
            @test dimension == 20
        end

        @testset "p₁ = x² + 1, p₂ = x² + x + 1, p₃ = x¹³ + 1" begin
            rel = [s * t * u + 1, s^2 + 1, t^2 + t + 1, u^13 + 1]

            relations = vcat(inv, rel)

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(Matrix, A, relations)

            # verify the relations are satisfied for low dimensions
            verify(relations, matrices)

            # test dimension
            @test dimension == 24
        end

        @testset "p₁ = x² + 1, p₂ = x² + x + 1, p₃ = x¹⁷ + 1" begin
            rel = [s * t * u + 1, s^2 + 1, t^2 + t + 1, u^17 + 1]

            relations = vcat(inv, rel)

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(Matrix, A, relations)

            # verify the relations are satisfied for low dimensions
            verify(relations, matrices)

            # test dimension
            @test dimension == 32
        end

        @testset "p₁ = x² + 1, p₂ = x³ + 1, p₃ = x³ + x + 1" begin
            rel = [s * t * u + 1, s^2 + 1, t^3 + 1, u^3 + u + 1]

            relations = vcat(inv, rel)

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(Matrix, A, relations)

            # verify the relations are satisfied for low dimensions
            verify(relations, matrices)

            # test dimension
            @test dimension == 9
        end

        @testset "p₁ = x² + 1, p₂ = x³ + 1, p₃ = x³ + x + 1" begin
            rel = [s * t * u + 1, s^2 + 1, t^3 + 1, u^3 + u + 1]

            relations = vcat(inv, rel)

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(Matrix, A, relations)

            # verify the relations are satisfied for low dimensions
            verify(relations, matrices)

            # test dimension
            @test dimension == 9
        end

        @testset "p₁ = x² + 1, p₂ = x³ + 1, p₃ = x⁵ + x² + 1" begin
            rel = [s * t * u + 1, s^2 + 1, t^3 + 1, u^5 + u^2 + 1]

            relations = vcat(inv, rel)

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(Matrix, A, relations)

            # verify the relations are satisfied for low dimensions
            verify(relations, matrices)

            # test dimension
            @test dimension == 25
        end
    end

    @testset "examples over GF(5)" begin
        𝔽₅ = GF(5)
        A, (s, t, u, s⁻¹, t⁻¹, u⁻¹) = free_associative_algebra(𝔽₅, ["s", "t", "u", "s⁻¹", "t⁻¹", "u⁻¹"])
        inv = [s * s⁻¹ - 1, s⁻¹ * s - 1, t * t⁻¹ - 1, t⁻¹ * t - 1, u * u⁻¹ - 1, u⁻¹ * u - 1]

        @testset "p₁ = x² - 1, p₂ = x² + x + 1, p₃ = x⁸ - 1" begin
            rel = [s * t * u - 1, s^2 - 1, t^2 + t + 1, u^8 - 1]

            relations = vcat(inv, rel)

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(Matrix, A, relations)

            # verify the relations are satisfied for low dimensions
            verify(relations, matrices)

            # test dimension
            @test dimension == 12
        end

        @testset "p₁ = x² - 1, p₂ = x² + x + 1, p₃ = x¹⁴ - 1" begin
            rel = [s * t * u - 1, s^2 - 1, t^2 + t + 1, u^14 - 1]

            relations = vcat(inv, rel)

            dimension = dimension_qa(A, relations)

            matrices = matrices_qa(Matrix, A, relations)

            # verify the relations are satisfied for low dimensions
            verify(relations, matrices)

            # test dimension
            @test dimension == 28
        end
    end

    M = 10
    @testset "families over GF(2) for m in 1:$M" begin
        𝔽₂ = GF(2)
        A, (s, t, u, s⁻¹, t⁻¹, u⁻¹) = free_associative_algebra(𝔽₂, ["s", "t", "u", "s⁻¹", "t⁻¹", "u⁻¹"])
        inv = [s * s⁻¹ - 1, s⁻¹ * s - 1, t * t⁻¹ - 1, t⁻¹ * t - 1, u * u⁻¹ - 1, u⁻¹ * u - 1]

        @testset "p₁ = x² - 1, p₂ = x² + x + 1, p₃ = xⁿ - 1" begin
            @testset "n = 2m" begin
                for n in 1:M
                    m = 2 * n
                    @testset "n = $m" begin
                        rel = [s * t * u + 1, s^2 + 1, t^2 + t + 1, u^m + 1]

                        relations = vcat(inv, rel)

                        dimension = dimension_qa(A, relations)
            
                        matrices = matrices_qa(Matrix, A, relations)
            
                        # verify the relations are satisfied for low dimensions
                        verify(relations, matrices)

                        # test dimension
                        @test dimension == 4 * n
                    end
                end
            end
            @testset "n = 2m + 1" begin
                for n in 1:M
                    m = 2 * n + 1
                    @testset "n = $m" begin
                        rel = [s * t * u + 1, s^2 + 1, t^2 + t + 1, u^m + 1]

                        relations = vcat(inv, rel)

                        dimension = dimension_qa(A, relations)
            
                        matrices = matrices_qa(Matrix, A, relations)
            
                        # verify the relations are satisfied for low dimensions
                        verify(relations, matrices)

                        # test dimension
                        @test dimension == 4 * n
                    end
                end
            end
        end
    end

    @testset "families over GF(5) for m in 1:$M" begin
        𝔽₅ = GF(5)
        A, (s, t, u, s⁻¹, t⁻¹, u⁻¹) = free_associative_algebra(𝔽₅, ["s", "t", "u", "s⁻¹", "t⁻¹", "u⁻¹"])
        inv = [s * s⁻¹ - 1, s⁻¹ * s - 1, t * t⁻¹ - 1, t⁻¹ * t - 1, u * u⁻¹ - 1, u⁻¹ * u - 1]
        @testset "p₁ = x² - 1, p₂ = x² + x + 1, p₃ = xⁿ - 1" begin
            @testset "n = 2m + 1" begin
                for n in 1:M
                    m = 2 * n + 1
                    @testset "n = $m" begin
                        rel = [s * t * u - 1, s^2 - 1, t^2 + t + 1, u^m - 1]

                        relations = vcat(inv, rel)

                        dimension = dimension_qa(A, relations)
            
                        matrices = matrices_qa(Matrix, A, relations)
            
                        # verify the relations are satisfied for low dimensions
                        verify(relations, matrices)

                        # test dimension
                        @test dimension == 0 broken = (n ∈ [1, 4, 7, 10])  # for these values the dimension is 2, even with the authors implementation
                    end
                end
            end
            @testset "n = 4m + 2" begin
                for n in 1:M
                    m = 4 * n + 2
                    @testset "n = $m" begin
                        rel = [s * t * u - 1, s^2 - 1, t^2 + t + 1, u^m - 1]

                        relations = vcat(inv, rel)

                        dimension = dimension_qa(A, relations)
            
                        matrices = matrices_qa(Matrix, A, relations)
            
                        # verify the relations are satisfied for low dimensions
                        verify(relations, matrices)

                        # test dimension
                        @test dimension == 8 * n + 4
                    end
                end
            end
            @testset "n = 4m" begin
                for n in 1:M
                    m = 4 * n
                    @testset "n = $m" begin
                        rel = [s * t * u - 1, s^2 - 1, t^2 + t + 1, u^m - 1]

                        relations = vcat(inv, rel)

                        dimension = dimension_qa(A, relations)
            
                        matrices = matrices_qa(Matrix, A, relations)
            
                        # verify the relations are satisfied for low dimensions
                        verify(relations, matrices)

                        # test dimension
                        @test dimension == 8 * n - 4
                    end
                end
            end
        end
    end
end
