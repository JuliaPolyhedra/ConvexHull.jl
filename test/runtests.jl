using ConvexHull, JuMP, Polyhedra, Test

@testset "Throws error in case of linearity" begin
    simplex = LiftedHRepresentation([0 1 0; 0 0 1; 1 -1 -1], BitSet([3]))
    @test_throws ErrorException double_description(simplex)
    linv = LiftedHRepresentation([1 1], BitSet([1]))
    @test_throws ErrorException double_description(linv)
end

include("polyhedron.jl")
