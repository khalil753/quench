# Abbreviations
using ForwardDiff

Fl, Vec, C = Float64, Vector, ComplexF64
VecCFl = Union{Vector{C}, Vector{Fl}}

struct NotImplmented <:Exception
  str::String
  function NotImplmented(str)
    println(str)    
  end
end

is_before_quench(X) = real(X[1]) <  0
is_after_quench(X)  = real(X[1]) >= 0

map_dict(f::Function, d::Dict)::Dict = Dict([(k, f(v)) for (k, v) in d])

get_crossed_derivative(f, ε) = (x,y) -> 1/4ε^2*(f(x + ε, y + ε) -  f(x - ε, y + ε) - f(x + ε, y - ε) + f(x - ε, y - ε))

function initialize_Ws(W_func, XA, XB)
  Ws = Dict()
  Ws["AB"] = DistributionWithTrajectories(W_func, XA, XB)
  Ws["AA"] = DistributionWithTrajectories(W_func, XA, XA)
  Ws["BB"] = DistributionWithTrajectories(W_func, XB, XB)
  return Ws
end

initialize_distributions(W_func, D_func, XA, XB) = (initialize_Ws(W_func, XA, XB), DistributionWithTrajectories(D_func, XA, XB))
function add_crossed_derivatives(Ws, D, ε_numeric_derivative)
  _gcd(f) = get_crossed_derivative(f, ε_numeric_derivative)
  return map_dict(_gcd, Ws), _gcd(D)
end

function negativity( ρ) 
  ρ22, ρ33 = real(ρ[2,2]), real(ρ[3,3])
  max(0, sqrt(abs2(ρ[1,4]) + ((ρ22 - ρ33)/2)^2) - (ρ22 + ρ33)/2)
end
function concurrence(ρ)
  ρ22, ρ33 = real(ρ[2,2]), real(ρ[3,3])
  # ρ22, ρ33 = max(real(ρ[2,2]), 0), max(real(ρ[3,3]), 0)
  2*max(0, abs(ρ[1,4]) - sqrt(ρ22*ρ33))
end