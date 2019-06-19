import Pkg

const POLYHEDRA_TEST = joinpath(dirname(dirname(pathof(Polyhedra))), "test")

include(joinpath(POLYHEDRA_TEST, "utils.jl"))
include(joinpath(POLYHEDRA_TEST, "solvers.jl")) # defines lp_solver
include(joinpath(POLYHEDRA_TEST, "polyhedra.jl"))

@testset "Polyhedra tests in $arith arithmetics" for arith in [:float, :exact]
    polyhedratest(ConvexHull.Library(arith, lp_solver), [
        "jumpsimplex", "nonfulldimensional", "simplex", "permutahedron",
        "board", "issue48", "empty"
    ])
end
