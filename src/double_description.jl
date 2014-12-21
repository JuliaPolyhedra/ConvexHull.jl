type DoubleDescription{T<:Real}
    A::Matrix{T}
    R::Vector{Vector{T}}
    K::Set{Int}
    adj::Vector{Vector{Bool}}
end
function DoubleDescription{T<:Real}(A::Matrix{T},R::Vector{Vector{T}},K::Set{Int})
    m = length(R)
    adj = [Array(Bool, i-1) for i in 1:m]
    d = rank(A)
    for i in 1:m
        for j in 1:(i-1)
            adj[i][j] = check_adjacency(A, R[i], R[j], d)
        end
    end
    return DoubleDescription(A,R,K,adj)
end

function initial_description{T<:Real}(A::Matrix{T})
    m,n = size(A)
    B = rref(A')
    # find pivots
    r = 1
    K = Set{Int}()
    for i in 1:m
        comp = zeros(n)
        comp[r] = 1
        if B[:,i] == comp
            r += 1
            push!(K, i)
        end
        r > n && break
    end
    Aₖ = A[collect(K),:]
    R = Aₖ \ eye(n,n)
    Rₖ = Vector{T}[R[:,i] for i in 1:n]
    return DoubleDescription(A,Rₖ,K)
end

function double_description{T<:Real}(A::Matrix{T})
    m, n = size(A)
    dd = initial_description(A)
    Kᶜ = setdiff(1:m, dd.K)
    while !isempty(Kᶜ)
        i = pop!(Kᶜ)
        Aᵢ = reshape(A[i,:], (n,))
        update!(dd, Aᵢ)
        push!(dd.K, i)
    end
    return dd
end

# use Lemma 8 from Fukuda (1996) to update the double description
function update!(dd::DoubleDescription, Aᵢ)
    for (j,k) in combinations(1:length(dd.R), 2)
        if isadjacent(dd,j,k)
            rⱼ, rₖ = dd.R[j], dd.R[k]
            tⱼ,tₖ = dot(Aᵢ,rⱼ), dot(Aᵢ,rₖ)
            if tⱼ > 0 && tₖ < 0
                addray!(dd, tⱼ*rₖ - tₖ*rⱼ)
            elseif tⱼ < 0 && tₖ > 0
                addray!(dd, tₖ*rⱼ - tⱼ*rₖ)
            end
        end
    end
    nothing
end

isadjacent(dd::DoubleDescription, j, k) = dd.adj[max(j,k)][min(j,k)]

function addray!(dd::DoubleDescription, r)
    m = length(dd.R)
    d = rank(dd.A)
    new_adj = Array(Bool, m)
    for i in 1:m
        new_adj[i] = check_adjacency(dd.A, r, dd.R[i], d)
    end
    push!(dd.adj, new_adj)
    push!(dd.R, r)
end

function check_adjacency(A, r, s, d)
    Z = active_sets(A, r, s)
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
