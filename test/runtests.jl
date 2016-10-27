using ConvexHull, JuMP, FactCheck

include("simplesquare.jl")
include("hypercube.jl")
include("simplex.jl")
include("ex1.jl")
include("infeasible.jl")
include("nonfulldimensional.jl")
include("crosspolytope.jl")

FactCheck.exitstatus()
