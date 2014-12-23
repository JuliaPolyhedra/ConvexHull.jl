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
    rankᴬ::Int
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
    Rₖ = [CountedVector{T}(R[:,i],i) for i in 1:n]
    dd = DoubleDescription(A,Rₖ,K,Dict{(Int,Int),Bool}(),n,rank(A))
    for i in 1:n
        for j in (i+1):n
            cache_adjacency!(dd, Rₖ[i], Rₖ[j])
        end
    end
    return dd
end

isredundant(x) = false

function double_description{T<:Real}(A::Matrix{T})
    A = [A; zeros(T,1,size(A,2))]
    A[end,1] = one(T)
    m, n = size(A)
    isredundant(A) && error("Cannot yet handle redundant systems")
    rank(A) < n && error("Cannot yet handle nontrivial lineality spaces")
    dd = initial_description(A)
    Kᶜ = setdiff(1:m, dd.K)
    while !isempty(Kᶜ)
        i = pop!(Kᶜ)
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
        val = minimum(findfirst(v .> ε))
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
            v = CountedVector(w,dd.num_rays)
            if sum(abs(w)) > n*ε && !is_approx_included(R⁰, vec(w)) && !is_approx_included(Rⁿᵉʷ, vec(w))
                dd.num_rays += 1
                push!(Rⁿᵉʷ, v)
            end
        end
    end
    dd.R = vcat(R⁺, R⁰, Rⁿᵉʷ)
    for r in dd.R, s in Rⁿᵉʷ
        r.id == s.id && continue
        cache_adjacency!(dd, r, s)
    end

    push!(dd.K, i)
    nothing
end

function partition_rays{T<:Real}(R::Vector{CountedVector{T}}, a::Vector{T})
    R⁺, R⁰, R⁻ = CountedVector{T}[], CountedVector{T}[], CountedVector{T}[]
    n = length(a)
    for r in R
        if sum(abs(vec(r))) < 10n*eps()
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

# isadjacent(dd, r, s) = check_adjacency(dd.A, r, s, rank(dd.A))
isadjacent(dd, r, s) = dd.adj[extrema([r.id,s.id])]

function cache_adjacency!(dd, r, s)
    d = dd.rankᴬ
    Z = active_sets(dd, r, s)
    if length(Z) < d - 2
        return dd.adj[extrema([r.id,s.id])] = false
    end
    dd.adj[extrema([r.id,s.id])] = (rank(dd.A[Z,:]) == d - 2)
end

function active_sets(dd, r, s)
    A = dd.A
    Ar = A*vec(r)
    Z = Int[]
    m = size(A,1)
    cnt = 0
    for i in 1:m
        if abs(Ar[i]) <= ε
            cnt += 1
        end
    end
    cnt < dd.rankᴬ - 2 && return Z # shortcurcuit early to avoid the secont mat-vec multiplication
    As = A*vec(s)
    sizehint!(Z,m)
    for i in 1:m
        if abs(Ar[i]) <= ε && abs(As[i]) <= ε
            push!(Z, i)
        end
    end
    return Z
end
