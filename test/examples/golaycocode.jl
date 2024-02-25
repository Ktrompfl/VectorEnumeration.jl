### group presentations for golay cocode from  
# S.A. Linton, Constructing matrix representations of finitely presented groups,
# Journal of Symbolic Computation, Volume 12, Issues 4–5, 1991,Pages 427-438, ISSN 0747-7171,
# https://doi.org/10.1016/S0747-7171(08)80095-8. (https://www.sciencedirect.com/science/article/pii/S0747717108800958)
𝔽₂ = GF(2)

@testset "M₂₃ - Golay cocode representation" begin
    A, (a, b, c, d, e, f) = free_associative_algebra(𝔽₂, ["a", "b", "c", "d", "e", "f"])
    M = free_module(A, 1)
    R = [a^2 - 1, b^2 - 1, c^2 - 1, d^2 - 1, e^2 - 1, f^2 - 1, 
        (a * b)^3 - 1, (a * e)^4 - 1, (b * c)^5 - 1, (c * d)^3 - 1, (c * e)^3 - 1, (c * f)^4 - 1,
        (d * f)^3 - 1, (e * f)^6 - 1, (a * c)^2 - 1, (a * d)^2 - 1, (a * f)^2 - 1, (b * d)^2 - 1, (b * e)^2 - 1, (b * f)^2 - 1,
        (d * e)^2 - 1, a * (c * f)^2 - 1, b * (e * f)^3 - 1, (e * a * b)^3 - 1, (b * c * e)^5 - 1, (a * e * c * d)^4 - 1, (b * c * e * f)^4 - 1
    ]
    W = [M([a - 1]), M([b - 1]), M([c - 1]), M([d - 1]), M([e - 1]), M([1 + f + f * e * a * b + f * d + f * d * c + f * d * c * e + f * d * c * e * f])]

    matrices = matrices_qm(Matrix, A, M, R, W)
    dimension = size(matrices[1], 1)

    verify(R, matrices)


    # test dimension
    @test dimension == 11
end

@testset "M₂₄ - Golay cocode representation" begin
    A, (a, b, c, d, e, f, g) = free_associative_algebra(𝔽₂, ["a", "b", "c", "d", "e", "f", "g"])
    M = free_module(A, 1)
    R = [ a^2 - 1, b^2 - 1, c^2 - 1, d^2 - 1, e^2 - 1, f^2 - 1, g^2 - 1,
        (a * b)^3 - 1, (b * c)^3 - 1, (b * f)^4 - 1, (c * d)^10 - 1, (d * g)^3 - 1,
        (d * e)^3 - 1, a * c * a * c - 1, a * d * a * d - 1, a * e * a * e - 1, a * f * a * f - 1, a * g * a * g - 1, b * d * b * d - 1, b * e * b * e - 1, b * g * b * g - 1, c * e * c * e - 1, c * f * c * f - 1,
        c * g * c * g - 1, d * f * d * f - 1, e * f * e * f - 1, e * g * e * g - 1, f * g * f * g - 1, a * (c * d)^5 - 1, a * (c * d * e)^5 - 1, f * (c * d * g)^5 - 1,
        e * (b * c * f)^3 - 1, (a * b * f)^3 - 1
    ]
    W = [M([a * b - 1]), M([a * c - 1]), M([d - 1]), M([g * e * a - 1]), M([e * d * f * b * c * d * e - 1]), M([1 + a + f + a * f + f * b + a * f * b + f * b * a + f * b * f])]

    matrices = matrices_qm(Matrix, A, M, R, W)
    dimension = size(matrices[1], 1)

    verify(R, matrices)

    # test dimension
    @test dimension == 12
end

@testset "2M₁₂ - Golay cocode representation" begin
    A, (d, x, y, z, t) = free_associative_algebra(𝔽₂, ["d", "x", "y", "z", "t"])
    M = free_module(A, 1)
    R = [ d^2 - 1, x^2 - 1, y^2 - 1, z^2 - 1, t^2 - 1,
        d * x * d * x - 1, d * y * d * y - 1, d * z * d * z - 1, d * t * d * t - 1, x * y * x * y - 1, x * z * x * z - 1, y * z * y * z - 1, (x * t)^3 - 1,
        (y * t)^5 - 1, (z * t)^5 - 1, d * (x * y * t)^5 - 1, d * (x * z * t)^5 - 1, (y * z * t)^5 - 1, d * (x * y * z * t)^5 - 1
    ]
    W = [M([x - 1]), M([t * y * t - 1]), M([t * z * t - 1]), M([t * y * z * y * t - 1]), M([t + y + z + y * z + t * x - 1])]

    dimension = dimension_qm(A, M, R, W)

    # test dimension
    @test dimension == 1  # according to the paper this should be 6, but even the authors implementation results in a 1 dimensional representation
end
