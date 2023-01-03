# WightmanFunction params
const space_time = "quench"

if     space_time=="quench"  const with_derivative_coupling = true
elseif space_time=="rindler" const with_derivative_coupling = true 
elseif space_time=="flat"    const with_derivative_coupling = false end

# Switching function params
const σ = 1.0
const switching_func_name = "cos4"
const switching_function_center_A = 0.5σ
const switching_function_center_B = 0.5σ

# Quench regulator
const b = 0.01

# Detector frequencies and Initial positions
nΩ = 20 # Number of frquencies/initial conditions to iterate over
nχ = 10
const χ0A = 0.5σ     
if     space_time=="quench"  Ωs, χ0Bs = LinRange(-1/σ, 30/σ, nΩ), LinRange(χ0A + 0.5σ, χ0A + 1.5σ, nχ)
elseif space_time=="rindler" Ωs, χ0Bs = LinRange(-1/σ, 30/σ, nΩ), LinRange(χ0A + 0.5σ, χ0A + 1.5σ, nχ)
elseif space_time=="flat"    Ωs, χ0Bs = LinRange(-3/σ,  3/σ, n) , LinRange(0.5σ      , 2σ        , n) end

# Coupling strength of detector
const λ = 1.0

# Delta for numerical derivation 
const ε_numeric_derivative = 1e-3

# Integrator params
if     space_time=="flat"    const initial_τs, final_τs = [switching_function_center_A - 1σ, -switching_function_center_B - 1σ], 
                                                          [switching_function_center_A + 1σ, switching_function_center_A + 1σ]
elseif space_time=="rindler" const initial_τs, final_τs = [-5σ, -5σ], [5σ, 5σ]
elseif space_time=="quench"  const initial_τs, final_τs = [max(σ*1e-2, switching_function_center_A - 5σ),
                                                          max(σ*1e-2, switching_function_center_B - 5σ)], 
                                                         [switching_function_center_A + 5σ, 
                                                          switching_function_center_B + 5σ] end
const int_tol = 1e-4
const maxevals = 220000

# Complex contour params
const ε_contour = 1e-2

params = Dict("Space Time"          => space_time,
              "Derivative Coupling" => with_derivative_coupling)
              ""
