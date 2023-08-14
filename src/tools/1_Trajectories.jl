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

@def_structequal struct AcceleratedTrajectory2D <: AbstractTrajectory
  χ0::Fl
end
(X::AcceleratedTrajectory2D)(τ) = [τ/abs(X.χ0), X.χ0]
