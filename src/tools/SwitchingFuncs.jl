_cos4( τ)::Complex = abs(real(τ)) <= 1/2 ? cos(π*τ)^4 : 0.0
_gauss(τ)::Complex = exp(-τ^2/2) 

switching_funcs = Dict("cos4"  => _cos4,
                       "gauss" => _gauss,
                       "one"   => x -> 1.0)

function initialize_switching_funcs(switching_func_name, σ, switching_function_center_A, switching_function_center_B)  
    χA(τ) = switching_funcs[switching_func_name]((τ - switching_function_center_A)/σ)
    χB(τ) = switching_funcs[switching_func_name]((τ - switching_function_center_B)/σ)
    [χA, χB]
end