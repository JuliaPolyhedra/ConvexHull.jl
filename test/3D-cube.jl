using CDD, JuMP

m = Model()
@defVar(m, 0 <= x[1:3] <= 1)

A, b = homogeneous_system(m)
dd = double_description([b A])

canonicalize!(dd)
