module CDD

#macro dd_ccall(func, args...)
#    f = "dd_$func"
#    quote
#        ccall(($f,libcdd), $(args...))
#    end
#end
#
#set_global_constants() = @dd_ccall(set_global_constants, ())
#set_global_constants()

export homogeneous_system, write_ine, read_ext

include("high_level.jl")
include("jump_interface.jl")

end
