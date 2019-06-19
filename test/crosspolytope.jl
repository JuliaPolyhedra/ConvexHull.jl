for N in 2:9
    facts("Cross polytope in $N dimensions") do
        m = Model()

        @variable(m, x[1:N])

        for k in 0:N
            for S in combinations(1:N, k)
                Sᶜ = setdiff(1:N, S)
                @constraint(m, sum(x[i] for i in S) - sum(x[i] for i in Sᶜ) ≤ 1)
            end
        end

        v, r = get_extrema(m)

        @fact length(v) --> 2N
        @fact length(r) --> 0

        for i in 1:N
            ex = zeros(N)
            ex[i] = 1
            @fact is_approx_included(v, ex) --> true
            ex[i] = -1
            @fact is_approx_included(v, ex) --> true
        end
    end
end
