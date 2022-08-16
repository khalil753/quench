_cos4(τ)  = abs(real(τ)) <= 1 ? cos(π*τ/2)^4 : 0
_gauss(τ) = exp(-τ^2/2σ) 

χs = Dict("cos4"  => _cos4,
          "gauss" => _gauss,
          "one"   => x -> 1.0)


