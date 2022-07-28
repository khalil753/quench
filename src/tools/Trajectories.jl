include("misc.jl")

abstract type AbstractTrajectory <: Function end

struct QuenchTrajectory <: AbstractTrajectory
  χ0::Fl
  y0::Fl
  b ::Fl 
  QuenchTrajectory(χ0::Fl, b::Fl) = new(χ0, asinh(χ0/b), b)
end
(X::QuenchTrajectory)(τ::Fl)::Vec{<:Fl} = τ > 0 ? [τ/X.χ0, X.y0] : [τ, X.χ0]

struct InertialTrajectory <: AbstractTrajectory
    x0::Fl
    y0::Fl
    z0::Fl
end
(X::InertialTrajectory)(τ::Fl)::Vec{<:Fl} = [τ, X.x0, X.y0, X.z0]

struct AcceleratedTrajectory <: AbstractTrajectory
  χ0::Fl
  y0::Fl
  z0::Fl
end
(X::AcceleratedTrajectory)(τ::Fl)::Vec{<:Fl} = [X.χ0*sinh(τ/X.χ0), X.χ0*cosh(τ/X.χ0), X.y0, X.z0]