module ConvexHull

using JuMP, GenericLinearAlgebra, RowEchelon, Combinatorics

using Polyhedra
using Printf: @printf

export homogeneous_system, read_ine, write_ine, read_ext, write_ext,
       double_description, canonicalize!,
       get_extrema, is_approx_included

const ε = 10eps()

include("double_description.jl")
include("jump_interface.jl")
include("readers_writers.jl")
include("polyhedron.jl")

end
