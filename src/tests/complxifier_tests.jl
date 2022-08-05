include("../quench.jl")

f(x) = 1/(x[1])
new_dist_func(x) = x[1]^2
f_c = 
g(x) = f_c([x[1], 0])
hcubature(g, [-5], [5], rtol=int_tol, maxevals=100000)[1]

