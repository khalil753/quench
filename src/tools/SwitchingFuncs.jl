_cos4(τ)::ComplexF64 = abs(real(τ)) <= 1 ? cos(π*τ/2)^4 : 0
_gauss(τ)::ComplexF64 = exp(-τ^2/2σ) 

χs = Dict("cos4"  => _cos4,
          "gauss" => _gauss,
          "one"   => x -> 1.0)


