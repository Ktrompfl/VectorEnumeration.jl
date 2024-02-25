using Documenter, Logging, Test
using AbstractAlgebra, SparseArrays, VectorEnumeration, LinearAlgebra

# DocMeta.setdocmeta!(VectorEnumeration, :DocTestSetup, :(using Nemo, VectorEnumeration); recursive=true)

include("util.jl")  # load common testing functions

@testset "VectorEnumeration.jl" begin
    # test docs
    doctest(VectorEnumeration; manual=true)

    @testset "Examples: VectorEnumeration" begin
        include("examples/permutation.jl")
        include("examples/golaycocode.jl")
        # these tests take very long to run
        # include("examples/cosets.jl")
        include("examples/polycyclics.jl")
    end

    @testset "Unittests: VectorEnumeration" begin
        include("unittests/algebratest.jl")
        include("unittests/relationtest.jl")
        # TODO: these tests need some clean up
        # include("unittests/processtest.jl")
        # include("unittests/vectortest.jl")
    end
end
