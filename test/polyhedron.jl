const polyhedra_test = joinpath(Pkg.dir("Polyhedra"), "test")

include(joinpath(polyhedra_test, "utils.jl"))
include(joinpath(polyhedra_test, "solvers.jl"))
include(joinpath(polyhedra_test, "polyhedra.jl"))

@testset "Polyhedra tests in $arith arithmetics" for arith in [:float, :exact]
    polyhedratest(ConvexHullLibrary(arith), ["jumpsimplex", "nonfulldimensional",
                                             "simplex", "permutahedron", "board", "issue48", "empty"])
end
