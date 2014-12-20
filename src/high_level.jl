function extrema{T<:Real}(A::Matrix{T}, b::Vector{T})
    
end

function facets{T<:Real}(vertices::Vector{Vector{T}}, rays::Vector{Vector{T}})
    
end

function write_ine(str::String, A::Matrix{Int}, b::Vector{Int})
    fp = open(str, "w")
    println(fp, "H-representation")
    println(fp, "begin")
    m, n = size(A)
    @assert n == length(b)
    println(fp, "$m $n integer")
    for i in 1:m
        println(fp, join([b[i], A[i,:]], " "))
    end
    println(fp, "end")
    close(fp)
    nothing
end

function write_ine(str::String, A::Matrix{Float64}, b::Vector{Float64})
    fp = open(str, "w")
    println(fp, "H-representation")
    println(fp, "begin")
    m, n = size(A)
    @assert m == length(b)
    println(fp, "$m $n integer")
    for i in 1:m
        @printf(fp, "%E ", b[i])
        for j in 1:n
            @printf(fp, "%E ", A[i,j])
        end
        println(fp)
    end
    println(fp, "end")
    close(fp)
    nothing
end
