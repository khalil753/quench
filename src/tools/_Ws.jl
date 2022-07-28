include("misc.jl")
include("Trajectories.jl")

struct NotImplmented <:Exception
  str::String
  function NotImplmented(str)
    println(str)    
  end
end

function W₀(z::T,z′::T)::C where {T<:Fl}
  iε = im*1e-1
  - im/8*sign(z - z′) - 1/(4π)*(log(cosh((z + z′)/2)) + log(sinh(abs(z - z′ + iε)/2)))
end

is_before_quench(X) = X[1] <  0
is_after_quench(X)  = X[1] >= 0

function _W_quench(X::Vec{T}, X′::Vec{T})::C where T<:Fl
  if is_after_quench(X) && is_after_quench(X′) 
    (η,y), (η′,y′) = X, X′
    return W₀(η - y, η′ - y′) + W₀(η + y, η′ + y′)
  elseif is_before_quench(X) && is_before_quench(X′)
    Δt, Δx = X - X′ 
    return -1/4π*(log(abs(Δt^2 - Δx^2))) - im/8*(sign(Δt+Δx) + sign(Δt-Δx))
  else 
    # return 0
    throw(NotImplmented("I haven't implememted the case where x and x′ belong to different spacetime patches")) 
  end
end

function _W_flat_spacetime(X::Vec{T}, X′::Vec{T})::C where T<:Fl
  ε = 1e-1
  Δt, Δx = X[1] - X′[1], X[2:end] - X′[2:end]
  -1/(4*π^2*((Δt - im*ε)^2 - sum(Δx.^2)))
end

_D(X::Vec{<:Fl}, X′::Vec{<:Fl})::C = if (X[1] > X′[1]) _W(X,X′) else _W(X′,X) end

struct CorrFuncWithTrajectories <: Function
  _G::Function
  X::AbstractTrajectory; X′::AbstractTrajectory 
end
CorrFuncWithTrajectories(_G, χ0::Fl, χ0′::Fl, b::Fl) = $correlation_func_name(_G, Trajectory(χ0,b), Trajectory(χ0′,b))

(G::CorrFuncWithTrajectories)(x::Vec{<:Fl}, x′::Vec{<:Fl})::C = G._G(x, x′)  
(G::CorrFuncWithTrajectories)(τ::Fl       , τ′::Fl       )::C = G(G.X(τ), W.X′(τ′))
