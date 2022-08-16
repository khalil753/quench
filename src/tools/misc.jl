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

function initialize_stuff()
  Wττ′s, Dττ′ = create_distributions(_Ws[space_time], _Ds[space_time], 
                                     χ0A, χ0B, b, 
                                     ε_numeric_derivative)
  χ(τ) = χs[switching_func_name]((τ)/σ)
  df = deform_funcs[deform_func_name]
  dist_func = distance_funcs[space_time]
  return Wττ′s, Dττ′, χ, df, dist_func
end


function get_Δτ(τs, deformation_function, pole_distance, ε)
  """This function gives me the complex part to add to τ, iΔτ"""
  d = pole_distance(τs)
  if d^2 <= ε^2  Δτ = ε*deformation_function(d/ε)
  else           Δτ = 0.0 end
end

function complexify(f, deformation_function, pole_distance, ε, get_∇Δτ) 
  function deformed_f(τs)::C
      Δτ          = get_Δτ(τs, deformation_function, pole_distance, ε)
      i∇Δτ, i∇Δτ′ = im.*get_∇Δτ(τs)     
      τ, τ′ = τs[1], τs[2] 
      # return f([τ - im*Δτ, τ′])*(1 - i∇Δτ)
      if Δτ > 0 return f([τ - im*Δτ, τ′])*(1 - i∇Δτ)
      else      return f([τ        , τ′])*(1 - i∇Δτ) end
  end
  return deformed_f
end

function complexify(f, deformation_function, pole_distance, ε) 
  _get_Δτ(τs) = get_Δτ(τs, deformation_function, pole_distance, ε)
  get_∇Δτ(τs) = ForwardDiff.gradient(_get_Δτ, τs)
  return complexify(f, deformation_function, pole_distance, ε, get_∇Δτ)
end

complexify(f, ε) = τs -> f([τs[1] - im*ε, τs[2] + im*ε])

# function get_crossed_derivative(f::Function)::Function
#   _f(xs) = f(xs[1], xs[2])
#   Rf, If = real∘_f, imag∘_f
#   fxy(x, y) = ForwardDiff.hessian(Rf, [x,y])[1,2] + im*ForwardDiff.hessian(If, [x,y])[1,2]
# end

get_crossed_derivative(f, ε) = (x,y) -> 1/4ε^2*(f(x + ε, y + ε) -  f(x - ε, y + ε) - f(x + ε, y - ε) + f(x - ε, y - ε))

function create_distributions(_W, _D, χ0A, χ0B, b, ε_numeric_derivative)
  XA, XB = QuenchTrajectory(χ0A, b), QuenchTrajectory(χ0B, b)
  Ws = Dict()
  Ws["AB"] = DistributionWithTrajectories(_W, XA, XB)
  Ws["AA"] = DistributionWithTrajectories(_W, XA, XA)
  Ws["BB"] = DistributionWithTrajectories(_W, XB, XB)

  _gcd(f) = get_crossed_derivative(f,ε_numeric_derivative)
  Wττ′s = map_dict(_gcd, Ws)
  Dττ′  = _gcd(DistributionWithTrajectories(_D, XA, XB))
  return  Wττ′s,  Dττ′
end

function negativity( ρ) 
  ρ22, ρ33 = real(ρ[2,2]), real(ρ[3,3])
  max(0, sqrt(abs2(ρ[1,4]) + ((ρ22 - ρ33)/2)^2) - (ρ22 + ρ33)/2)
end
function concurrence(ρ)
  ρ22, ρ33 = real(ρ[2,2]), real(ρ[3,3])
  2*max(0, abs(ρ[1,4]) - sqrt(ρ22*ρ33))
end
