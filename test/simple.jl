using JuMP, CDD

m = Model()
@defVar(m, 0 <= x[i=1:3] <= i/2)
for i in 1:3
    @addConstraint(m, dot(rand(3),x) <= i)
end

write_ine(m, "simple.ine")
