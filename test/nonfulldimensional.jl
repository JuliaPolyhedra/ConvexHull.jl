facts("Non full-dimensional") do
    m = Model()
    @variable(m, x[1:3] ≥ 1)
    @constraints(m, begin
        x[1] == 2
        x[2] ≤ 2
    end)

    v, r = get_extrema(m)

    @fact is_approx_included(v, [2,2,1]) --> true
    @fact is_approx_included(v, [2,1,1]) --> true
    @fact is_approx_included(r, [0,0,1]) --> true
end
