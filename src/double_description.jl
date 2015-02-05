type DoubleDescription{T<:Real}
    A::Matrix{T}
    R::Vector
    K::Set{Int}
    L #Vector{LinkedList}
end

type CountedVector{T<:Real}
    v::Vector{T}
    Av::Vector{T}
    dd::DoubleDescription{T}

    function CountedVector(v::Vector{T},dd::DoubleDescription{T})
        canonicalize!(v)
        new(v,dd.A*v,dd)
    end
end
CountedVector{T<:Real}(v::Vector{T},dd::DoubleDescription{T}) = CountedVector{T}(v,dd)

Base.show(io::IO, v::CountedVector) = show(io, v.v)

Base.vec(c::CountedVector) = c.v

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
    cK  = sort(collect(K))
    cKᶜ = setdiff(1:m,cK) 
    A = vcat(A[cK,:],
             A[cKᶜ,:])
    Aₖ = A[1:n,:]
    d = rank(Aₖ)
    R = Aₖ \ eye(n,n)
    dd = DoubleDescription(A,
                           CountedVector{T}[],
                           Set{Int}(1:n),
                           [(CountedVector{T},CountedVector{T})[] for i in 1:m])
    Rₖ = [CountedVector{T}(R[:,i],dd) for i in 1:n]
    dd.R = Rₖ
    for i in 1:n, j in 1:n
        i == j && continue
        r, s = Rₖ[i], Rₖ[j]
        if satisfies_cₖ(r,s)
            if is_adjacent(dd, r.Av[1:n], s.Av[1:n], d)
                conditional_edge_store!(dd, r, s, mii(r), n)
            end
        end
    end
    return dd
end

function mii(r)
    val = findfirst(r.Av) do x
        x < -ε
    end
    val == 0 && (val = length(r.Av)+1)
    return val
end

function satisfies_cₖ(r,s)
    (k = mii(r)) < mii(s) || return false
    abs(s.Av[k]) ≤ ε && return false
    return true
end

function is_adjacent(dd, Ar, As, d)
    Z = active_sets(dd, Ar, As)
    if length(Z) < d - 2
        return false
    elseif length(intersect(Z,dd.K)) ≥ d - 2
        return true
    else
        return (rank(dd.A[sort(collect(dd.K)),:][Z,:]) == d - 2)
    end
end

function conditional_edge_store!(dd, r, s, k, i)
    Z = active_sets(dd, r.Av, s.Av)
    incl = any((i+1):(k-1)) do j
        j in Z
    end
    if !incl
        # if is_adjacent(dd, r, s)
            # dd.L[k] = list(dd.L[k], (r,s))
            push!(dd.L[k], (r,s))
        # end
    end
end

function double_description{T<:Real}(A::Matrix{T})
    A = [zeros(T,1,size(A,2)); A]
    A[1,1] = one(T)
    m, n = size(A)
    dd = initial_description(A)
    nK = length(dd.K)
    for i in (nK+1):m
        println("iteration $i of $m")
        # length(dd.K) == 31 && break
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
        tmp = filter(v) do x
            abs(x) > ε
        end
        val = abs(minimum(tmp))
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
    for (r,s) in dd.L[i]
        if r ∈ R⁺
            @assert s ∈ R⁻
            w = dot(Aᵢ,vec(r))*vec(s) - dot(Aᵢ,vec(s))*vec(r)
        elseif r ∈ R⁻
            @assert s ∈ R⁺
            w = dot(Aᵢ,vec(s))*vec(r) - dot(Aᵢ,vec(r))*vec(s)
        else
            error("why is this happening??")
        end
        v = CountedVector(w,dd)
        # if sum(abs(w)) > n*ε && 
           # !is_approx_included(R⁰,   vec(w)) && 
           # !is_approx_included(Rⁿᵉʷ, vec(w))
            push!(Rⁿᵉʷ, v)
        # end
    end
    dd.R = vcat(R⁺, R⁰, Rⁿᵉʷ)
    push!(dd.K, i)
    d = rank(dd.A[sort(collect(dd.K)),:])
    for r in vcat(R⁰, Rⁿᵉʷ), s in vcat(R⁰, Rⁿᵉʷ)
        if satisfies_cₖ(r,s)
            conditional_edge_store!(dd, r, s, mii(r), i)
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
