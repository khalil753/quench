using StructEquality
include("misc.jl")

abstract type AbstractTrajectory <: Function end
@def_structequal struct QuenchTrajectory <: AbstractTrajectory
  χ0::Fl
  y0::Fl
  b ::Fl
  sqrt_chi_b::Fl
  QuenchTrajectory(χ0::Fl, b::Fl) = new(χ0, asinh(χ0/b), b, sqrt(χ0^2 + b^2))
end
(X::QuenchTrajectory)(τ) = real(τ) > 0 ? [τ/X.sqrt_chi_b, X.y0] : [τ, X.χ0]

@def_structequal struct InertialTrajectory <: AbstractTrajectory
    x0::Fl
    y0::Fl
    z0::Fl
end
(X::InertialTrajectory)(τ) = [τ, X.x0, X.y0, X.z0]
  
@def_structequal struct AcceleratedTrajectory <: AbstractTrajectory
  χ0::Fl
  x0::Fl
  y0::Fl
  z0::Fl
end
(X::AcceleratedTrajectory)(τ) = [X.χ0*sinh(τ/abs(X.χ0)), X.χ0*cosh(τ/abs(X.χ0)) + X.x0 - X.χ0, X.y0, X.z0]

get_γ(X::AbstractTrajectory)    = (println("You haven't created the gamma function for this trajectory"), throw(Exception))
get_γ(X::InertialTrajectory)    = τ -> 1
get_γ(X::QuenchTrajectory)      = τ -> 1/X.χ0
# get_γ(X::QuenchTrajectory)      = τ -> cosh(τ/X.χ0)
# get_γ(X::AcceleratedTrajectory) = τ -> cosh(τ/X.χ0)
# get_γ(X::AcceleratedTrajectory) = τ -> 1

lorentz_distance(X, X′) = (X[1] - X′[1])^2 - sum((X[2:end] - X′[2:end]).^2)

function quench_distance(X, X′) 
  if is_after_quench(X) && is_after_quench(X′) || is_before_quench(X) && is_before_quench(X′)
    lorentz_distance(X, X′)
  else 
    warn("I haven't implemented the distance between spacetime points when they are 
                            \nin different patches of the quench")
    0.0
  end
end

distance_funcs = Dict("flat"   => lorentz_distance,
                      "quench" => quench_distance)