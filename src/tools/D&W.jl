include("misc.jl")
include("Trajectories.jl")

struct NotImplmented <:Exception
  str::String
  function NotImplmented(str)
    println(str)    
  end
end

function W₀(z, z′, ε)
  iε = im*ε
  - 1/(4π)*(log(cosh((z + z′)/2)) + log(sinh(z - z′ - iε)/2)) # + log(sinh(abs(z - z′ - iε)/2)))
end

function _W_quench(X::Vec{C}, X′::Vec{C}, ε::Fl)
  if is_after_quench(X) && is_after_quench(X′) 
    (η,y), (η′,y′) = X, X′
    return W₀(η - y, η′ - y′, ε) + W₀(η + y, η′ + y′, ε)
  elseif is_before_quench(X) && is_before_quench(X′)
    Δt, Δx = X - X′ 
    iε = im*ε
    return -1/4π*(log(abs((Δt-iε)^2 - Δx^2))) - im/8*(sign(Δt+Δx) + sign(Δt-Δx))
  else 
    # return 0
    throw(NotImplmented("I haven't implememted the case where x and x′ belong to different spacetime patches")) 
  end
end

function _W_flat_spacetime(X::Vec{C}, X′::Vec{C}, ε::C)
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
  _G::Function
  X ::AbstractTrajectory
  X′::AbstractTrajectory 
  ε ::Float64
  Δf::Function
end
function DistributionWithTrajectories(_G::Function, χ0::Fl, χ0′::Fl, b::Fl, ε ::Float64, Δf::Function) 
  """ Initialization for quench trajectories"""
  DistributionWithTrajectories(_G, QuenchTrajectory(χ0,b), QuenchTrajectory(χ0′,b), ε, Δf)
end

function (G::DistributionWithTrajectories)(X::Vec{C}, X′::Vec{C})::C
  d = abs(distance(X, X′))
  ε = G.ε
  d <= ε && G._G(X, X′, ε*G.Δf(d/ε))  
  else      return G._G(X, X′) end

(G::DistributionWithTrajectories)(τ::Fl, τ′::Fl)::C = G(G.X(τ), G.X′(τ′))

