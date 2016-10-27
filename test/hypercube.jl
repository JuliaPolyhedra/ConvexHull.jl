for N in 2:11
    facts("Hypercube in $N dimensions") do
        m = Model()
        @variable(m, 0 ≤ x[1:N] ≤ 1)

        v, r = get_extrema(m)

        @fact length(v) --> 2^N
        @fact length(r) --> 0

        for k in 0:N
            for p in combinations(1:N, k)
                ex = zeros(N)
                ex[p] = ones(length(p))
                @fact ex in v --> true
            end
        end
    end
end
