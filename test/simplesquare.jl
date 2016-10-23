facts("A simple square") do
    m = Model()
    @variable(m, 0 ≤ x[1:2] ≤ 1)

    v, r = get_extrema(m)

    @fact length(v) --> 4
    @fact length(r) --> 0
    @fact [0,0] in v --> true
    @fact [1,0] in v --> true
    @fact [0,1] in v --> true
    @fact [1,1] in v --> true
end
