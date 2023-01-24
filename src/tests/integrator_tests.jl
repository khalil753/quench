using HCubature

# In this test the integral should be iπ and it is
d=1000
iε = im*1e-8
f(t) = 1/(t - iε)
f_vec(t) = f(t[1])
hcubature(f_vec, [-d], [d], maxevals=100000)# , rtol=rtol)