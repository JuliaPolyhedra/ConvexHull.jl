type CountedVector{T<:Real}
    v::Vector{T}
    id::Int

    CountedVector(v::Vector{T},id::Int) = new(canonicalize!(v),id)
end
CountedVector{T<:Real}(v::Vector{T},id::Int) = CountedVector{T}(v,id)

Base.vec(c::CountedVector) = c.v

type DoubleDescription{T<:Real}
    A::Matrix{T}
    R::Vector{CountedVector{T}}
    K::Set{Int}
    adj::Dict{(Int,Int),Bool}
    num_rays::Int
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
    Aₖ = A[sort(collect(K)),:]
    R = Aₖ \ eye(n,n)
    Rₖ = [CountedVector{T}(R[:,i],i) for i in 1:n]
    dd = DoubleDescription(A,Rₖ,K,Dict{(Int,Int),Bool}(),n)
    for i in 1:n
        Ar = Aₖ*vec(Rₖ[i])
        for j in (i+1):n
            As = Aₖ*vec(Rₖ[j])
            id = extrema([Rₖ[i].id,Rₖ[j].id])
            cache_adjacency!(dd, Aₖ, n, Ar, As, id)
        end
    end
    return dd
end

function double_description{T<:Real}(A::Matrix{T})
    A = [zeros(T,1,size(A,2)); A]
    A[1,1] = one(T)
    m, n = size(A)
    dd = initial_description(A)
    Kᶜ = setdiff(1:m, dd.K)
    while !isempty(Kᶜ)
        i = pop!(Kᶜ)
        println("iteration $(length(dd.K)+1) of $m")
        update!(dd, i)
    end
    return dd
end

function is_approx_included(haystack, needle)
    n = length(needle)
    diff = zeros(n)
    for h in haystack
        diff = norm(vec(h)-needle)
        diff < n*ε && return true
    end
    return false
end

function canonicalize!{T<:Real}(v::Vector{T})
    n = length(v)
    val = abs(v[1])
    if val < ε
        val = abs(v[findfirst(abs(v) .> ε)])
    end
    for i in 1:n
        v[i] = v[i] / val
    end
    v
end

# use Lemma 8 from Fukuda (1996) to update the double description
function update!{T<:Real}(dd::DoubleDescription{T}, i)
    m, n = size(dd.A)
    Aᵢ = reshape(dd.A[i,:], (n,))
    Rⁿᵉʷ = CountedVector{T}[]
    R⁺, R⁰, R⁻ = partition_rays(dd.R, Aᵢ)
    for r in R⁺, s in R⁻
        if isadjacent(dd,r,s)
            w = dot(Aᵢ,vec(r))*vec(s) - dot(Aᵢ,vec(s))*vec(r)
            v = CountedVector(w,dd.num_rays+1)
            if sum(abs(w)) > n*ε && 
               !is_approx_included(R⁰,   vec(w)) && 
               !is_approx_included(Rⁿᵉʷ, vec(w))
                dd.num_rays += 1
                push!(Rⁿᵉʷ, v)
            end
        end
    end
    dd.R = vcat(R⁺, R⁰, Rⁿᵉʷ)
    push!(dd.K, i)
    Aₖ = dd.A[sort(collect(dd.K)),:]
    # should really add a test right about here to ensure
    # that old rays do not become adjacent...I think this 
    # can only happen if both v,w ∈ R⁰
    d = rank(Aₖ)
    for s in Rⁿᵉʷ
        As = Aₖ*vec(s)
        for r in dd.R
            r.id == s.id && continue
            Ar = Aₖ*vec(r)
            id = extrema([r.id, s.id])
            cache_adjacency!(dd, Aₖ, d, Ar, As, id)
        end
    end
    nothing
end

function partition_rays{T<:Real}(R::Vector{CountedVector{T}}, a::Vector{T})
    R⁺, R⁰, R⁻ = CountedVector{T}[], CountedVector{T}[], CountedVector{T}[]
    n = length(a)
    for r in R
        if sum(abs(vec(r))) < n*ε
            println("have a zero vector!")
            continue
        end
        val = dot(a,vec(r))
        if val > ε
            push!(R⁺, r)
        elseif val < -ε
            push!(R⁻, r)
        else
            push!(R⁰, r)
        end
    end
    return R⁺, R⁰, R⁻
end

isadjacent(dd, r, s) = dd.adj[extrema([r.id,s.id])]

function cache_adjacency!(dd, Aₖ, d, Ar, As, id)
    Z = active_sets(dd, Ar, As)
    if length(Z) < d - 2
        return (dd.adj[id] = false)
    end
    if length(intersect(Z,dd.K)) ≥ d - 2
        return (dd.adj[id] = true)
    end
    dd.adj[id] = (rank(Aₖ[Z,:]) == d - 2)
end

function active_sets(dd, Ar, As)
    Z = Int[]
    m = length(Ar)
    sizehint!(Z,m)
    for i in 1:m
        if abs(Ar[i]) ≤ ε && abs(As[i]) ≤ ε
            push!(Z, i)
        end
    end
    return Z
end
