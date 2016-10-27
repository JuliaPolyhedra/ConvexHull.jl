facts("The ex1 example") do
    m = Model()
    @variable(m, x)
    @variable(m, y)
    @constraints(m, begin
        12 + 2x - y ≥ 0
        -6 - x + 2y ≥ 0
        -3 + x + y ≥ 0
        1 + x ≥ 0
    end)

    v, r = get_extrema(m)

    @fact length(r) --> 2
    @fact length(v) --> 3

    for ex in Vector[[ 0, 3],
        [-1, 4],
        [-1,10]]
        @fact is_approx_included(v, ex) --> true
    end

    for ex in Vector[[1,2],
        [2,1]]
        @fact is_approx_included(r, ex) --> true
    end
end
