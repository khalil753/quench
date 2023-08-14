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
    warn("Watch out, I haven't implemented the case where x and x′ belong to different spacetime patches")
    return 0
  end
end

function _W_flat_2d(X, X′)
  Δt, Δx = X[1] - X′[1], X[2:end] - X′[2:end]
  -1/(4*π)*log(Δt^2 - sum(Δx.^2))
end

_W_rindler2 = _W_flat_2d

function _W_rindler(X, X′)
  (η, χ), (η′, χ′) = X, X′
  return -1/(4*π)*log(χ^2 + χ′^2 - 2χ*χ′*cosh(η - η′))
end

function _time_order(W::Function)
  """This function creates a new time_ordered_W which is time oredered (duh)"""
  time_ordered_W(X, X′) = real(X[1]) > real(X′[1]) ? W(X,X′) : W(conj(X′),conj(X)) 
  return time_ordered_W
end

_Ws = make_Ws_dict()
_Ds = map_dict(_time_order, _Ws)

@def_structequal struct DistributionWithTrajectories <: Function
  _G::Function
  X ::Function
  X′::Function
end
DwT = DistributionWithTrajectories

(G::DwT)(X::Vec   , X′::Vec   ) = G._G(  X   ,   X′    ) 
(G::DwT)(τ::Number, τ′::Number) = G._G(G.X(τ), G.X′(τ′)) 