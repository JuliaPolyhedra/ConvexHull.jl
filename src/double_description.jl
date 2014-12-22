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
    Aₖ = A[sort(collect(K)),:]
    R = Aₖ \ eye(n,n)
    Rₖ = Vector{T}[R[:,i] for i in 1:n]
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

# use Lemma 8 from Fukuda (1996) to update the double description
function update!{T<:Real}(dd::DoubleDescription{T}, i)
    m, n = size(dd.A)
    Aᵢ = reshape(dd.A[i,:], (n,))
    Rⁿᵉʷ = Vector{T}[]
    R⁺, R⁰, R⁻ = partition_rays(dd.R, Aᵢ)
    for r in R⁺, s in R⁻
        if isadjacent(dd,r,s)
            v = dot(Aᵢ,r)*s - dot(Aᵢ,s)*r
            if sum(abs(v)) > 10n*eps()
                push!(Rⁿᵉʷ, v)
            else
                error("Somehow we added a zero vector?")
            end
        end
    end
    # println("Rⁿᵉʷ = $(Rⁿᵉʷ)")
    # println("R⁺   = $(R⁺)")
    # println("R⁰   = $(R⁰)")
    # println("R⁻   = $(R⁻)")
    dd.R = vcat(Rⁿᵉʷ, R⁺, R⁰)
    push!(dd.K, i)
    # canonicalize!(dd)
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
        if val > 0#10eps()
            push!(R⁺, r)
        elseif val < 0#-10n*eps()
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
