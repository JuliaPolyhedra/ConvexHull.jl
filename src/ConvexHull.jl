module ConvexHull

using JuMP, LinearAlgebra, RowEchelon, Combinatorics

export homogeneous_system, read_ine, write_ine, read_ext, write_ext,
       double_description, canonicalize!,
       get_extrema, is_approx_included

using Polyhedra

const ε = 10eps()

include("double_description.jl")
include("jump_interface.jl")
include("readers_writers.jl")
include("polyhedron.jl")

end
