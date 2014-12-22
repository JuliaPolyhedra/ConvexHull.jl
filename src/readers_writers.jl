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

function read_ext(fname::String)
    fp = open(fname, "r")
    vertices = Vector[]
    rays = Vector[]

    while true
        line = strip(readline(fp), '\n')
        line == "begin" && break
    end
    line = strip(readline(fp), '\n')
    (m,n,typ) = filter(x->!isempty(x), split(line, " "))
    line = strip(readline(fp), '\n')
    while line != "end"
        res  = filter(x->!isempty(x), split(line, " "))
        println("res = $res")
        nums = map(float, res)
        elem = try convert(Vector{Int}, nums) catch nums end
        if elem[1] == 1 # vertex
            push!(vertices, elem[2:end])
        else
            push!(rays, elem[2:end])
        end
        line = strip(readline(fp), '\n')
    end
    close(fp)
    return vertices, rays
end
