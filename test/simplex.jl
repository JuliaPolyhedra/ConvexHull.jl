using CDD, JuMP

m = Model()
@defVar(m, x[1:6] >= 0)
@addConstraint(m, sum{x[i], i=1:6} <= 1)

A, b = homogeneous_system(m)
dd = double_description([b A])

canonicalize!(dd)
