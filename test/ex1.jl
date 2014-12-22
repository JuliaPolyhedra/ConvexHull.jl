using CDD, JuMP

m = Model()
@defVar(m, x)
@defVar(m, y)
@addConstraints(m, begin
    12 + 2x - y ≥ 0
    -6 - x + 2y ≥ 0
    -3 + x + y ≥ 0
    1 + x ≥ 0
end)

A, b = homogeneous_system(m)
dd = double_description([b A])

CDD.canonicalize!(dd)
