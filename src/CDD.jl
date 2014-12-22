module ConvexHull

using JuMP

export homogeneous_system, write_ine, read_ext, double_description, canonicalize!, 
       extrema, is_approx_included

const ε = 10eps()

include("double_description.jl")
include("JuMP_interface.jl")
include("readers_writers.jl")

end
