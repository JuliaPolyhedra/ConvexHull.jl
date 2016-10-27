# The 6D simplex
for N in 1:3:51
    facts("Simplex in $N dimensions") do
        N = 6
        m = Model()
        @variable(m, x[1:N] >= 0)
        @constraint(m, sum{x[i], i=1:N} == 1)

        v, r = get_extrema(m)

        @fact length(v) --> N
        @fact length(r) --> 0

        for k in 1:N
            ex = zeros(N)
            ex[k] = 1
            @fact ex in v --> true
        end
    end
end

# The 6D simplex, plus the origin
for N in 1:3:39
    facts("Simplex with the origin in $N dimensions") do
        m = Model()
        @variable(m, x[1:N] >= 0)
        @constraint(m, sum{x[i], i=1:N} ≤ 1)

        v, r = get_extrema(m)

        @fact length(v) --> N+1
        @fact length(r) --> 0

        @fact zeros(N) in v --> true
        for k in 1:N
            ex = zeros(N)
            ex[k] = 1
            @fact ex in v --> true
        end
    end
end
