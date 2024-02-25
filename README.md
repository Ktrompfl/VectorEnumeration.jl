# VectorEnumeration
VectorEnumeration is a Julia package for constructing matrix representations of finitely generated algebras over fields using vector enumeration.

---

## Example of Usage

```Julia

julia> using Nemo

julia> using VectorEnumeration

julia> field = QQ
Rational field

julia> A, (a, b, c) = free_associative_algebra(QQ, ["a", "b", "c"])
(Free associative algebra on 3 indeterminates over QQ, AbstractAlgebra.Generic.FreeAssAlgElem{QQFieldElem}[a, b, c])

julia> inv = [a, b, c]  # all generators are involutions
3-element Vector{AbstractAlgebra.Generic.FreeAssAlgElem{QQFieldElem}}:
 a
 b
 c

julia> rel = [(a*b)^3, (b*c)^3, (a*c)^2]  # fixing relators
3-element Vector{AbstractAlgebra.Generic.FreeAssAlgElem{QQFieldElem}}:
 a*b*a*b*a*b
 b*c*b*c*b*c
 a*c*a*c

julia> rank = 1
1

julia> sub = [(a*b*c - 1, )]
1-element Vector{Tuple{AbstractAlgebra.Generic.FreeAssAlgElem{QQFieldElem}}}:
 (a*b*c - 1,)

julia> (generator, matrices), (dimension, ntotal), (images, preimages) = vector_enumeration(A, inv, rel, rank, sub)
((AbstractAlgebra.Generic.FreeAssAlgElem{QQFieldElem}[a, b, c], SparseArrays.SparseMatrixCSC{QQFieldElem, Int64}[sparse([2, 1, 4, 3, 6, 5], [1, 2, 3, 4, 5, 6], QQFieldElem[1, 1, 1, 1, 1, 1], 6, 6), sparse([5, 3, 2, 6, 1, 4], [1, 2, 3, 4, 5, 6], QQFieldElem[1, 1, 1, 1, 1, 1], 6, 6), sparse([3, 4, 1, 2, 6, 5], [1, 2, 3, 4, 5, 6], QQFieldElem[1, 1, 1, 1, 1, 1], 6, 6)]), (6, 7), (nothing, SparseArrays.SparseVector{QQFieldElem, Int64}[  [1]  =  1]))

julia> matrices[1]  #a
6×6 SparseArrays.SparseMatrixCSC{QQFieldElem, Int64} with 6 stored entries:
 ⋅  1  ⋅  ⋅  ⋅  ⋅
 1  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  1  ⋅  ⋅
 ⋅  ⋅  1  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  1
 ⋅  ⋅  ⋅  ⋅  1  ⋅

julia> matrices[2]  #b
6×6 SparseArrays.SparseMatrixCSC{QQFieldElem, Int64} with 6 stored entries:
 ⋅  ⋅  ⋅  ⋅  1  ⋅
 ⋅  ⋅  1  ⋅  ⋅  ⋅
 ⋅  1  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  1
 1  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  1  ⋅  ⋅

julia> matrices[3]  #c
6×6 SparseArrays.SparseMatrixCSC{QQFieldElem, Int64} with 6 stored entries:
 ⋅  ⋅  1  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  1  ⋅  ⋅
 1  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  1  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  1
 ⋅  ⋅  ⋅  ⋅  1  ⋅

julia> dimension
6

julia> ntotal
7
```

## References
[S. Linton, Vector Enumeration Programs, version 3 (1993)]("https://github.com/gap-packages/ve").

[S. Linton, On Vector Enumeration, (1993)]("https://www.sciencedirect.com/science/article/pii/002437959390245J")
