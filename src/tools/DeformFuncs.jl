_cos(τ)  = abs(τ) <= 0.5 ? cos(π*τ)   : 0.0
_cos2(τ) = abs(τ) <= 0.5 ? cos(π*τ)^2 : 0.0
# _cos4(τ) = abs(τ) <= 0.5 ? cos(π*τ)^4 : 0.0

function _triangle(t) 
    if     -1 <  t <= 0 return (1 + t)
    elseif  0 <= t <= 1 return (1 - t) 
    else                        0 
    end
end

# _∇cos4(τ) = abs(τ) <  1 ? π/2*cos(π*τ/2)^3*sin(π*τ/2) : 0

deform_funcs  = Dict("cos4"     => _cos4,
                     "cos2"     => _cos2,
                     "cos"      => _cos,
                     "triangle" => _triangle)

# ∇deform_funcs = Dict("cos4"     => _∇cos4)
