for N in 2:8
    facts("Infeasible in $N dimensions") do
        m = Model()
        @variable(m, 0 ≤ x[1:N] ≤ 1)
        @constraint(m, x[1] ≥ 2)

        r, v = get_extrema(m)

        @fact isempty(r) --> true
        @fact isempty(v) --> true
    end
end
