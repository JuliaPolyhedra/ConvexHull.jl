using JuMP

function homogeneous_system(model::JuMP.Model)
    A = full(JuMP.prepConstrMatrix(model))
    c, lb, ub = JuMP.prepProblemBounds(model)
    l, u = model.colLower, model.colUpper

    m, n = size(A)
    @assert m == length(lb) == length(ub)
    @assert model.nlpdata == nothing
    @assert isempty(model.quadconstr)
    @assert isempty(model.sosconstr)

    C = Array(Float64, 0, n)
    b = Float64[]
    for i in 1:m
        if !isinf(lb[i])
            C = [C; A[i,:]']
            push!(b, -lb[i])
        end
        if !isinf(ub[i])
            C = [C; -A[i,:]']
            push!(b, ub[i])
        end
    end

    F = eye(n, n)
    B = Array(Float64, 0, n)
    d = Float64[]
    for i in 1:n
        if !isinf(l[i])
            B = [B; F[i,:]']
            push!(d, -l[i])
        end
        if !isinf(u[i])
            B = [B; -F[i,:]']
            push!(d, u[i])
        end
    end

    return [C;B], [b;d]
end

function write_ine(model::JuMP.Model, fname::String)
    A, b = homogeneous_system(model)
    try
        write_ine(fname, convert(Matrix{Int},A), convert(Vector{Int},b))
    catch InexactError
        try
            write_ine(fname, convert(Matrix{Int},2A), convert(Vector{Int},2b))
        catch InexactError
            try
                write_ine(fname, convert(Matrix{Int},4A), convert(Vector{Int},4b))
            catch InexactError
                write_ine(fname, A, b)
            end
        end
    end
end

function get_extrema(model::JuMP.Model)
    A, b = homogeneous_system(model)
    dd = double_description([b A])

    vertices, rays = Vector{Float64}[], Vector{Float64}[]
    for cv in dd.R
        r = vec(cv)
        # if abs(r[1]) > 10eps()
        if abs(r[1] - 1) < ε
            push!(vertices, r[2:end])
            # push!(vertices, r[2:end] ./ r[1])
        else
            @assert abs(r[1]) < ε
            push!(rays, r[2:end] ./ minimum(r[find(r)]))
        end
    end
    return vertices, rays
end

