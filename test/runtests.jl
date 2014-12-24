using ConvexHull, JuMP, Base.Test

# A simple square
let
    m = Model()
    @defVar(m, 0 ≤ x[1:2] ≤ 1)

    v, r = get_extrema(m)

    @test length(v) == 4
    @test length(r) == 0
    @test [0,0] in v
    @test [1,0] in v
    @test [0,1] in v
    @test [1,1] in v
end

# the hypercube
for N in 2:11
    let
        m = Model()
        @defVar(m, 0 ≤ x[1:N] ≤ 1)

        v, r = get_extrema(m)

        @test length(v) == 2^N
        @test length(r) == 0
        
        for k in 0:N
            for p in combinations(1:N, k)
                ex = zeros(N)
                ex[p] = ones(length(p))
                @test ex in v
            end
        end
    end
end

# The 6D simplex
for N in 1:3:51
    let
        N = 6
        m = Model()
        @defVar(m, x[1:N] >= 0)
        @addConstraint(m, sum{x[i], i=1:N} == 1)

        v, r = get_extrema(m)

        @test length(v) == N
        @test length(r) == 0
        
        for k in 1:N
            ex = zeros(N)
            ex[k] = 1
            @test ex in v
        end
    end
end

# The 6D simplex, plus the origin
for N in 1:3:39
    let
        m = Model()
        @defVar(m, x[1:N] >= 0)
        @addConstraint(m, sum{x[i], i=1:N} ≤ 1)

        v, r = get_extrema(m)

        @test length(v) == N+1
        @test length(r) == 0
        
        @test zeros(N) in v
        for k in 1:N
            ex = zeros(N)
            ex[k] = 1
            @test ex in v
        end
    end
end

# The ex1 example
let
    m = Model()
    @defVar(m, x)
    @defVar(m, y)
    @addConstraints(m, begin
        12 + 2x - y ≥ 0
        -6 - x + 2y ≥ 0
        -3 + x + y ≥ 0
        1 + x ≥ 0
    end)

    v, r = get_extrema(m)

    @test length(r) == 2
    @test length(v) == 3

    for ex in Vector[[ 0, 3],
                     [-1, 4],
                     [-1,10]]
        @test is_approx_included(v, ex)
    end

    for ex in Vector[[1,2],
                     [2,1]]
        @test is_approx_included(r, ex)
    end
end

# infeasible
for N in 2:8
    let
        m = Model()
        @defVar(m, 0 ≤ x[1:N] ≤ 1)
        @addConstraint(m, x[1] ≥ 2)

        r, v = get_extrema(m)

        @test isempty(r)
        @test isempty(v)
    end
end

# non full-dimensional 
let
    m = Model()
    @defVar(m, x[1:3] ≥ 1)
    @addConstraints(m, begin
        x[1] == 2
        x[2] ≤ 2
    end)

    v, r = get_extrema(m)

    @test is_approx_included(v, [2,2,1])
    @test is_approx_included(v, [2,1,1])
    @test is_approx_included(r, [0,0,1])
end

# cross polytope
for N in 2:8
    let
        m = Model()

        @defVar(m, x[1:N])

        for k in 0:N
            for S in combinations(1:N, k)
                Sᶜ = setdiff(1:N, S)
                @addConstraint(m, sum{x[i], i in S} - sum{x[i], i in Sᶜ} ≤ 1)
            end
        end

        v, r = get_extrema(m)

        @test length(v) == 2N
        @test length(r) == 0

        for i in 1:N
            ex = zeros(N)
            ex[i] = 1
            @test is_approx_included(v, ex)
            ex[i] = -1
            @test is_approx_included(v, ex)
        end
    end
end
