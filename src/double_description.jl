type DoubleDescription{T<:Real}
    A::Matrix{T}
    R::Vector{Vector{T}}
    K::Set{Int}
    adj::Dict{(Int,Int),Bool}
end
function DoubleDescription{T<:Real}(A::Matrix{T},R::Vector{Vector{T}},K::Set{Int})
    m = length(R)
    adj = Dict{(Int,Int),Bool}()
    d = rank(A)
    for i in 1:m
        for j in 1:(i-1)
            adj[(i,j)] = check_adjacency(A, R[i], R[j], d)
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

function double_description{T<:Real}(A::Matrix{T})
    m, n = size(A)
    if rank(A) ≥ n
        dd = initial_description(A)
    else
        error("Cannot yet handle nontrivial lineality spaces")
        # dd = smaller_initial_description(A)
    end
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
function update!{T<:Real}(dd::DoubleDescription{T}, Aᵢ::Vector{T})
    Rⁿᵉʷ = Vector{T}[]
    R⁺, R⁰, R⁻ = partition_rays(dd.R, Aᵢ)
    n = length(Aᵢ)
    for r in R⁺, s in R⁻
        if isadjacent(dd,r,s)
            v = dot(Aᵢ,r)*s - dot(Aᵢ,s)*r
            if abs(sum(v)) > 10n*eps()
                push!(Rⁿᵉʷ, v)
            end
        end
    end
    # println("Rⁿᵉʷ = $(Rⁿᵉʷ)")
    # println("R⁺   = $(R⁺)")
    # println("R⁰   = $(R⁰)")
    # println("R⁻   = $(R⁻)")
    dd.R = vcat(Rⁿᵉʷ, R⁺, R⁰)
    nothing
end

function partition_rays{T<:Real}(R::Vector{Vector{T}}, a::Vector{T})
    R⁺, R⁰, R⁻ = Vector{T}[], Vector{T}[], Vector{T}[]
    n = length(a)
    for r in R
        val = dot(a,r)
        if abs(sum(r)) < 10n*eps()
            println("have a zero vector!")
            continue
        end
        if val > 10eps()
            push!(R⁺, r)
        elseif val < -10n*eps()
            push!(R⁻, r)
        else
            push!(R⁰, r)
        end
    end
    return R⁺, R⁰, R⁻
end

# isadjacent(dd::DoubleDescription, j, k) = dd.adj[max(j,k)][min(j,k)]
isadjacent(dd, j, k) = true

function addray!(dd::DoubleDescription, r)
    m = length(dd.R)
    d = rank(dd.A)
    # for i in 1:m
    #     new_adj[i] = check_adjacency(dd.A, r, dd.R[i], d)
    # end
    # push!(dd.adj, new_adj)
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

function canonicalize!(dd::DoubleDescription)
    m, n = size(dd.A)
    for i in 1:length(dd.R)
        r₁ = dd.R[i][1]
        if r₁ != 0
            dd.R[i] ./= r₁
        end
    end
end

function purge_invalid_rays{T<:Real}(dd::DoubleDescription{T})
    Aₖ = dd.A[collect(dd.K),:]
    R = Vector{T}[]
    sizehint!(R, length(dd.R))
    for r in dd.R
        s = Aₖ*r
        nonneg = all(s) do x
            x ≥ 0
        end
        nonneg && push!(R, r)
    end
    dd.R = R
    nothing
end
