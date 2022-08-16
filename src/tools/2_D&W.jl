using ForwardDiff
include("1_Trajectories.jl")

W₀(z::C , z′::C ) = - 1/4π*(log((cosh((z + z′))/2)*sinh(   (z - z′)/2))) 
W₀(z    , z′    ) = - 1/4π*(log((cosh((z + z′))/2)*sinh(abs(z - z′)/2))) - im/8*sign(z - z′)

function _W_quench(X, X′)
  if is_after_quench(X) && is_after_quench(X′) 
    (η,y), (η′,y′) = X, X′
    return W₀(η - y, η′ - y′) + W₀(η + y, η′ + y′)
  elseif is_before_quench(X) && is_before_quench(X′)
    Δt, Δx = X - X′ 
    return -1/4π*(log((Δt)^2 - Δx^2))
  else 
    warn("Watch out, I haven't implememted the case where x and x′ belong to different spacetime patches")
    return 0
  end
end

function _W_flat_spacetime(X, X′)
  Δt, Δx = X[1] - X′[1], X[2:end] - X′[2:end]
  -1/(4*π^2*(Δt^2 - sum(Δx.^2)))
end

function _time_order(W::Function)
  """This function creates a new time_ordered_W which is time oredered (duh)"""
  time_ordered_W(X, X′) = real(X[1]) > real(X′[1]) ? W(X,X′) : W(X′,X) 
  return time_ordered_W
end

_Ws = Dict("quench" => _W_quench,
           "flat"   => _W_flat_spacetime)
_Ds = map_dict(_time_order, _Ws)

struct DistributionWithTrajectories <: Function
  _G::Function
  X ::AbstractTrajectory
  X′::AbstractTrajectory 
end

(G::DistributionWithTrajectories)(X::Vec   , X′::Vec   ) = G._G(  X   ,   X′    ) 
(G::DistributionWithTrajectories)(τ::Number, τ′::Number) = G._G(G.X(τ), G.X′(τ′)) 