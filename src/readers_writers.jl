function write_ine(str::String, A::Matrix{Int}, b::Vector{Int})
    fp = open(str, "w")
    println(fp, "H-representation")
    println(fp, "begin")
    m, n = size(A)
    @assert m == length(b)
    println(fp, "$m $(n+1) integer")
    for i in 1:m
        println(fp, join([b[i] A[i,:]], " "))
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
    println(fp, "$m $(n+1) real")
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

function write_ext(str::String, ext::Vector)
    fp = open(str, "w")
    println(fp, "V-representation")
    println(fp, "begin")
    m = length(ext)
    n = length(ext[1])
    println(fp, "$m $n real")
    for i in 1:m
        vv = ext[i]
        @assert length(vv) == n
        println(fp, join(vv, " "))
    end
    println(fp, "end")
    println(fp, "hull")
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
    elem = []
    while line != "end"
        res  = filter(x->!isempty(x), split(line, " "))
        nums = map(float, res)
        try
            elem = convert(Vector{Int}, nums)
        catch
            elem = nums
        end
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

function read_ine(fname::String)
    fp = open(fname, "r")

    while true
        line = strip(readline(fp), '\n')
        line == "begin" && break
    end
    line = strip(readline(fp), '\n')
    (mm,nn,typ) = filter(x->!isempty(x), split(line, " "))
    m,n = int(mm), int(nn)
    line = strip(readline(fp), '\n')
    A = Array(Float64,0,n-1)
    b = Float64[]
    while line != "end"
        res  = filter(x->!isempty(x), split(line, " "))
        elem = map(float, res)
        push!(b, elem[1])
        A = [A; -elem[2:end]']
        line = strip(readline(fp), '\n')
    end
    @assert size(A) == (m,n-1)
    close(fp)
    return b, A
end
