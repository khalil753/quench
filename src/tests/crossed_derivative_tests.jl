using ForwardDiff
include("../quench.jl")

function gcd(f::Function)::Function
    _f(xs) = f(xs[1], xs[2])
    Rf, If = real∘_f, imag∘_f
    fxy(x, y) = ForwardDiff.hessian(Rf, [x,y])[1,2] + im*ForwardDiff.hessian(If, [x,y])[1,2]
end
  
function gcd2(f::Function, ε::Fl)::Function
    fxy(x, y) = (f(x+ε, y+ε) - f(x+ε, y-ε) - f(x-ε, y+ε) + f(x-ε, y-ε))/(4*ε^2)
end
  
XA, XB = QuenchTrajectory(χ0A, b), QuenchTrajectory(χ0B, b)
W = DistributionWithTrajectories(_Ws["quench"] , XA, XB)
A, B = gcd(W), gcd2(W, 7.5e-7)

println(A(1,2))
println(B(1,2))