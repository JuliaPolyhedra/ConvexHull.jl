module ConvexHull

using JuMP

export homogeneous_system, read_ine, write_ine, read_ext, write_ext,
       double_description, canonicalize!,
       get_extrema, is_approx_included

importall Base

const ε = 10eps()

include("double_description.jl")
include("JuMP_interface.jl")
include("readers_writers.jl")

end
