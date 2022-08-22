_cos4( τ)::Complex = abs(real(τ)) <= 1 ? cos(π*τ/2)^4 : 0.0
_gauss(τ)::Complex = exp(-τ^2/2) 

χs = Dict("cos4"  => _cos4,
          "gauss" => _gauss,
          "one"   => x -> 1.0)


