import Pkg

const POLYHEDRA_TEST = joinpath(dirname(dirname(pathof(Polyhedra))), "test")

include(joinpath(POLYHEDRA_TEST, "utils.jl"))
include(joinpath(POLYHEDRA_TEST, "config.jl"))
include(joinpath(POLYHEDRA_TEST, "solvers.jl")) # defines lp_solver
include(joinpath(POLYHEDRA_TEST, "jump.jl"))

include(joinpath(POLYHEDRA_TEST, "basic.jl"))
include(joinpath(POLYHEDRA_TEST, "simplex.jl"))
include(joinpath(POLYHEDRA_TEST, "permutahedron.jl"))
include(joinpath(POLYHEDRA_TEST, "board.jl"))
include(joinpath(POLYHEDRA_TEST, "docexample.jl"))
include(joinpath(POLYHEDRA_TEST, "issue48.jl"))
include(joinpath(POLYHEDRA_TEST, "empty.jl"))
include(joinpath(POLYHEDRA_TEST, "sparse.jl"))
include(joinpath(POLYHEDRA_TEST, "sparserect.jl"))
include(joinpath(POLYHEDRA_TEST, "support_function.jl"))
const misctests = Dict("basic" => basictest,
                       "doc" => doctest,
                       "simplex" => simplextest,
                       "permutahedron" => permutahedrontest,
                       "board" => boardtest,
                       "issue48" => issue48test,
                       "empty" => emptytest,
                       "sparse" => sparsetest,
                       "sparserect" => sparserecttest,
                       "support_function" => support_function_test)

@polytestset misc

include(joinpath(POLYHEDRA_TEST, "decompose.jl"))

polyhedratests = Dict("jump" => jumptest,
                      "misc" => misctest,
                      "decompose" => decomposetest)

@polytestset polyhedra true

@testset "Polyhedra tests in $arith arithmetics" for arith in [:float, :exact]
    polyhedratest(ConvexHull.Library(arith, lp_solver), [
        "jumpsimplex", "nonfulldimensional", "simplex", "permutahedron",
        "board", "issue48", "empty"
    ])
end
