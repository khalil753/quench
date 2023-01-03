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

function initialize_switching_funcs(switching_func_name, σ, switching_function_center_A, switching_function_center_B)  
  χA(τ) = switching_funcs[switching_func_name]((τ-switching_function_center_A)/σ)
  χB(τ) = switching_funcs[switching_func_name]((τ-switching_function_center_B)/σ)
  [χA, χB]
end

get_crossed_derivative(f, ε) = (x,y) -> 1/4ε^2*(f(x + ε, y + ε) -  f(x - ε, y + ε) - f(x + ε, y - ε) + f(x - ε, y - ε))

function initialize_Ws(W_func, XA, XB)
  Ws = Dict()
  Ws["AB"] = DistributionWithTrajectories(W_func, XB, XA)
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
  if ρ22 <= -1e-3 println("PB is quite negative so there could be a numerical problem: PB = $ρ22") end
  if ρ33 <= -1e-3 println("PA is quite negative so there could be a numerical problem: PA = $ρ33") end
  ρ22, ρ33 = max(ρ22, 0), max(ρ33, 0)
  max(0, sqrt(abs2(ρ[1,4]) + ((ρ22 - ρ33)/2)^2) - (ρ22 + ρ33)/2)
end
function concurrence(ρ)
  ρ22, ρ33 = real(ρ[2,2]), real(ρ[3,3])
  if ρ22 <= -1e-3 println("PB is quite negative so there could be a numerical problem: PB = $ρ22") end
  if ρ33 <= -1e-3 println("PA is quite negative so there could be a numerical problem: PA = $ρ33") end
  ρ22, ρ33 = max(ρ22, 0), max(ρ33, 0)
  2*max(0, abs(ρ[1,4]) - sqrt(ρ22*ρ33))
end