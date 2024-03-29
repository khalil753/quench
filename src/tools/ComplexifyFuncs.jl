# function get_Δτ(τs, deformation_function, pole_distance, ε)
#     """This function gives me the complex part to add to τ, iΔτ"""
#     d = pole_distance(τs)
#     if d^2 <= ε^2  Δτ = ε*deformation_function(d/ε)
#     else           Δτ = 0.0 end
# end
  
# function complexify(f, deformation_function, pole_distance, ε, get_∇Δτ) 
#     function deformed_f(τs)::C
#         Δτ      = get_Δτ(τs, deformation_function, pole_distance, ε)
#         i∇Δτ, _ = im.*get_∇Δτ(τs)     
#         τ, τ′ = τs[1], τs[2] 
#         # return f([τ - im*Δτ, τ′])*(1 - i∇Δτ)
#         if Δτ > 0 return f([τ - im*Δτ, τ′])*(1 - i∇Δτ)
#         else      return f([τ        , τ′])*(1 - i∇Δτ) end
#     end
#     return deformed_f
# end

# function complexify(f, deformation_function, pole_distance, ε) 
#     _get_Δτ(τs) = get_Δτ(τs, deformation_function, pole_distance, ε)
#     get_∇Δτ(τs) = ForwardDiff.gradient(_get_Δτ, τs)
#     return complexify(f, deformation_function, pole_distance, ε, get_∇Δτ)
# end

complexify(f, ε) = (τ1,τ2) -> f(τ1 - im*ε, τ2 + im*ε)

function deform(f, ε) 
    iε = im*ε
    function new_f(X, Y)
        deformed_X = [X[1] - iε, X[2:end]...]
        deformed_Y = [Y[1] + iε, Y[2:end]...]
        f(deformed_X, deformed_Y)
    end
end