_cos4( τ)::Complex = abs(real(τ)) <= 1/2 ? cos(π*τ)^4 : 0.0
_gauss(τ)::Complex = exp(-τ^2/2) 

switching_funcs = Dict("cos4"  => _cos4,
                       "gauss" => _gauss,
                       "one"   => x -> 1.0)

function initialize_switching_funcs(switching_func_name, σ, ηcA_or_τcA, ηcB_or_τcB, using_ηs, χ0A, χ0B)
    if using_ηs 
         sqrtA, sqrtB = sqrt(χ0A^2 + b^2), sqrt(χ0B^2 + b^2)
         τcA, τcB = ηcA_or_τcA*sqrtA, ηcB_or_τcB*sqrtB 
    else τcA, τcB = ηcA_or_τcA      , ηcB_or_τcB       end
    χA(τ) = switching_funcs[switching_func_name]((τ - τcA)/σ)
    χB(τ) = switching_funcs[switching_func_name]((τ - τcB)/σ)
    [χA, χB]
end