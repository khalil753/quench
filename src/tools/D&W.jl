include("misc.jl")
include("Trajectories.jl")

struct NotImplmented <:Exception
  str::String
  function NotImplmented(str)
    println(str)    
  end
end

function W₀(z::T,z′::T, ε::Fl)::C where {T<:Fl}
  iε = im*ε
  - im/8*sign(z - z′) - 1/(4π)*(log(cosh((z + z′)/2)) + log(sinh(z - z′ - iε)/2)) # + log(sinh(abs(z - z′ - iε)/2)))
end

function _W_quench(X::Vec{T}, X′::Vec{T}, ε::Fl)::C where T::C
  iε = im*ε
  if is_after_quench(X) && is_after_quench(X′) 
    (η,y), (η′,y′) = X, X′
    return W₀(η - y, η′ - y′, ε) + W₀(η + y, η′ + y′, ε)
  elseif is_before_quench(X) && is_before_quench(X′)
    Δt, Δx = X - X′ 
    return -1/4π*(log(abs((Δt-iε)^2 - Δx^2))) - im/8*(sign(Δt+Δx) + sign(Δt-Δx))
  else 
    # return 0
    throw(NotImplmented("I haven't implememted the case where x and x′ belong to different spacetime patches")) 
  end
end

function _W_flat_spacetime(X::Vec{T}, X′::Vec{T}, ε::T)::C where T<:Fl
  iε = im*ε
  Δt, Δx = X[1] - X′[1], X[2:end] - X′[2:end]
  -1/(4*π^2*((Δt - iε)^2 - sum(Δx.^2)))
end

function time_order_dist(distrib_func::Function)
  function time_ordered_distrib_func(X::Vec{T}, X′::Vec{T}, ε::T)::C where T<: Fl
    X[1] > X′[1] ? distrib_func(X,X′,ε) : distrib_func(X′,X,ε) 
  end
  return time_ordered_distrib_func
end

_Ws = Dict("quench" => _W_quench,
           "flat"   => _W_flat_spacetime)
_Ds = map_dict(time_order_dist, _Ws)

struct DistributionWithTrajectories <: Function
  distrib_func::Function
  X ::AbstractTrajectory
  X′::AbstractTrajectory 
end
function DistributionWithTrajectories(distrib_func::Function, χ0::Fl, χ0′::Fl, b::Fl) 
  """ Initialization for quench trajectories"""
  DistributionWithTrajectories(distrib_func, QuenchTrajectory(χ0,b), QuenchTrajectory(χ0′,b))
end

(G::DistributionWithTrajectories)(X::Vec{<:Fl}, X′::Vec{<:Fl}, ε::Fl=0.0)::C = G.distrib_func(  X   ,   X′    , ε)  
(G::DistributionWithTrajectories)(τ::Fl       , τ′::Fl       , ε::Fl=0.0)::C = G.distrib_func(G.X(τ), G.X′(τ′), ε)

