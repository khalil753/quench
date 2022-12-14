using StructEquality
include("misc.jl")

abstract type AbstractTrajectory <: Function end
@def_structequal struct QuenchTrajectory <: AbstractTrajectory
  χ0::Fl
  y0::Fl
  b ::Fl
  QuenchTrajectory(χ0::Fl, b::Fl) = new(χ0, asinh(χ0/b), b)
end
(X::QuenchTrajectory)(τ) = real(τ) > 0 ? [τ/X.χ0, X.y0] : [τ, X.χ0]

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
(X::AcceleratedTrajectory)(τ) = [X.χ0*sinh(τ/X.χ0), X.χ0*cosh(τ/X.χ0) + X.x0 - X.χ0, X.y0, X.z0]

get_γ(X::AbstractTrajectory) = τ -> ((X(τ + 1e-5) - X(τ - 1e-5))/2e-5)[1]
get_γ(X::InertialTrajectory) = τ -> 1
get_γ(X::QuenchTrajectory)   = τ -> 1/X.χ0

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