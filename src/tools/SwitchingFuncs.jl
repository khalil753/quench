include("params.jl")

_cos(  τ::AbstractFloat)::Complex = abs(τ) <= σ ? cos(π*τ/(2σ)) : 0
_gauss(τ::AbstractFloat)::Complex = exp(-τ^2/(2σ)) 

χs = Dict("cos"   => _cos,
          "gauss" => _gauss,
          "one"   => x -> 1.0)


