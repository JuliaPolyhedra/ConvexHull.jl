const polyhedra_test = joinpath(Pkg.dir("Polyhedra"), "test")

include(joinpath(polyhedra_test, "alltests.jl"))

tests = Tuple{String, Function}[]
push!(tests, ("Hypercube in 2 dimensions", lib->hypercubetest(lib, 2)))
#push!(tests, ("Simplex in 2 dimensions", lib->simplextest(lib, 2))) # linearity
push!(tests, ("Simplex with the origin in 2 dimensions", lib->simplexorigtest(lib, 2)))
push!(tests, ("Cross Polytope in 2 dimensions", lib->crosspolytopetest(lib, 2)))
push!(tests, ("The ex1 example", ex1test))
push!(tests, ("Infeasible in 2 dimensions", lib->infeasibletest(lib, 2)))
#push!(tests, ("Non full-dimensional", nonfulldimensionaltest)) #linearity

@testset "Polyhedra tests in $arith arithmetics" for arith in [:float, :exact]
    runtests(ConvexHullLib(arith), tests)
end
