_cos4(τ)  = abs(τ) <= 1 ?     cos(π*τ/2)^4 : 0
_∇cos4(τ) = abs(τ) <  1 ? π/2*cos(π*τ/2)^3*sin(π*τ/2) : 0
deform_funcs  = Dict("cos4" => _cos4)

∇deform_funcs = Dict("cos4" => _∇cos4)