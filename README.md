# VectorEnumeration.jl

A [Julia](https://julialang.org) package for constructing matrix representations of finitely presented algebras over fields using vector enumeration.

--- 

## Introduction
VectorEnumeration.jl is an implementation of the vector enumeration algorithm, as presented by S.A. Linton in [[1]](https://github.com/Ktrompfl/VectorEnumeration.jl#references) and [[2]](https://github.com/Ktrompfl/VectorEnumeration.jl#references), and is based on the [latest available version](https://github.com/gap-packages/ve) of the authors implementation from 1996.

The package builds upon the types from [AbstractAlgebra.jl](https://github.com/Nemocas/AbstractAlgebra.jl) and provides additional functionality when [Hecke.jl](https://github.com/thofma/Hecke.jl) or [Oscar.jl](https://github.com/oscar-system/Oscar.jl) are loaded.

## Installation
VectorEnumeration.jl is no registered package yet and therefore must be installed from sources.
 
## Quick start
Here is an example of using VectorEnumeration.jl to compute the permutation matrix representation over $\mathbb{Q}$ and a $\mathbb{Q}$-base of the dihedral group of order 6 for the presentation $\langle a, b\ |\ a^3 = b^2 = (ab)^2 =1 \rangle$:

```jldoctest
julia> using AbstractAlgebra, SparseArrays, VectorEnumeration

julia> A, (x, y) = free_associative_algebra(QQ, [:x, :y]);

julia> R = [x^3 - 1, y^2 - 1, (x*y)^2 - 1];

julia> X, Y = matrices_qa(SparseMatrixCSC, A, R);

julia> X
6×6 SparseMatrixCSC{Rational{BigInt}, Int64} with 6 stored entries:
 ⋅  1  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  1  ⋅  ⋅  ⋅
 1  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  1
 ⋅  ⋅  ⋅  1  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  1  ⋅

julia> Y
6×6 SparseMatrixCSC{Rational{BigInt}, Int64} with 6 stored entries:
 ⋅  ⋅  ⋅  1  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  1  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  1
 1  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  1  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  1  ⋅  ⋅  ⋅

julia> base_qa(A, R)
6-element Vector{Vector{AbstractAlgebra.Generic.FreeAssAlgElem{Rational{BigInt}}}}:
 [1]
 [x]
 [x^2]
 [y]
 [x*y]
 [x^2*y]
```

## License
VectorEnumeration.jl is licensed under the MIT license; see [LICENSE](https://github.com/Ktrompfl/VectorEnumeration.jl/blob/main/LICENSE) for the full license text.

## [References]

[1] S.A. Linton, *Constructing matrix representations of finitely presented groups*,
Journal of Symbolic Computation, Volume 12, Issues 4–5, 1991, Pages 427-438, ISSN 0747-7171,
<https://doi.org/10.1016/S0747-7171(08)80095-8>.

[2] S.A. Linton, *On vector enumeration*, 
Linear Algebra and its Applications, Volume 192, 1993, Pages 235-248, ISSN 0024-3795,
<https://doi.org/10.1016/0024-3795(93)90245-J>.