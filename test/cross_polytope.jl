using CDD, JuMP

m = Model()

N = 3

@defVar(m, -1 <= x[1:N] <= 1)

for k in 0:N
    for S in combinations(1:N, k)
        Sᶜ = setdiff(1:N, S)
        @addConstraint(m, sum{x[i], i in S} - sum{x[i], i in Sᶜ} ≤ 1)
    end
end

A, b = homogeneous_system(m)
dd = double_description([b A])

CDD.canonicalize!(dd)
