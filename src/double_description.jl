type DoubleDescription{T<:Real}
    A::Matrix{T}
    R::Vector{Vector{T}}
    K::Set{Int}
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
    Rₖ = Vector{T}[canonicalize!(R[:,i]) for i in 1:n]
    return DoubleDescription(A,Rₖ,K)
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
        diff = norm(h-needle)
        diff < n*ε && return true
    end
    return false
end

function canonicalize!(v)
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
    Rⁿᵉʷ = Vector{T}[]
    R⁺, R⁰, R⁻ = partition_rays(dd.R, Aᵢ)
    for r in R⁺, s in R⁻
        if isadjacent(dd,r,s)
            v = dot(Aᵢ,r)*s - dot(Aᵢ,s)*r
            canonicalize!(v)
            if sum(abs(v)) > n*ε && !is_approx_included(R⁰, v) && !is_approx_included(Rⁿᵉʷ, v)
                push!(Rⁿᵉʷ, v)
            # else
                # error("Somehow we added a zero vector?")
            end
        end
    end
    dd.R = vcat(Rⁿᵉʷ, R⁺, R⁰)
    push!(dd.K, i)
    nothing
end

function partition_rays{T<:Real}(R::Vector{Vector{T}}, a::Vector{T})
    R⁺, R⁰, R⁻ = Vector{T}[], Vector{T}[], Vector{T}[]
    n = length(a)
    for r in R
        if sum(abs(r)) < 10n*eps()
            println("have a zero vector!")
            continue
        end
        val = dot(a,r)
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

isadjacent(dd, r, s) = check_adjacency(dd.A, r, s, rank(dd.A))


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
    for i in 1:m
        if abs(Ar[i]) <= ε && abs(As[i]) <= ε
            push!(Z, i)
        end
    end
    return Z
end
