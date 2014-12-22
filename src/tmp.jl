function smaller_initial_description{T<:Real}(A::Matrix{T})
    K = Set{Int}(1)
    Aₖ = A[1,:]
    m, n = size(A)
    R = Array(Vector{T}, n)
    R[1] = reshape(Aₖ, (n,))
    x = Array(Vector{T}, n-1)
    k = findfirst(Aₖ)
    for i in 1:(k-1)
        x[i] = zeros(T,n)
        x[i][i] = one(T)
    end
    for i in (k+1):n
        x[i-1] = zeros(T,n)
        x[i-1][i] = one(T)
    end
    for i in 1:(n-1)
        r = x[i]
        for j in 1:i
            r -= dot(R[j],x[i]) / dot(R[j],R[j]) * R[j]
        end
        R[i+1] = r
    end
    append!(R, -R[2:end])
    return DoubleDescription(A,R,K)
end


function check_adjacency(A, r, s, d)
    Z = active_sets(A, r, s)
    length(Z) < d - 2 && return false
    return rank(A[Z,:]) == d - 2
end

function active_sets(A, r, s)
    Ar = A*r
    As = A*s
    m = size(A,1)
    Z = Int[]
    sizehint!(Z,m)
    ϵ = 10eps()
    for i in 1:m
        if abs(Ar[i]) <= ϵ && abs(As[i]) <= ϵ
            push!(Z, i)
        end
    end
    return Z
end
