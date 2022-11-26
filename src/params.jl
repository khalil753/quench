# WightmanFunction params
const space_time = "flat"
if     space_time=="quench" const with_derivative_coupling = true
elseif space_time=="flat"   const with_derivative_coupling = false end

# Switching function params
const σ = 1.0
const switching_func_name = "gauss"
const switching_function_center = 0σ
# const switching_function_center_A = 0σ
# const switching_function_center_B = 0σ

# Quench regulator
const b = 0.01

# Detector frequencies and Initial positions
const χ0A = 1.0
if     space_time=="quench" Ωs, χ0Bs = LinRange(-1/σ, 15/σ, 8), LinRange(χ0A + 0.5σ, χ0A + 2.5σ, 8)
elseif space_time=="flat"   Ωs, χ0Bs = LinRange(-3/σ,  3/σ, 8), LinRange(0.5σ      , 2σ        , 8) end

# Coupling strength of detector
const λ = 1.0

# Delta for numerical derivation 
const ε_numeric_derivative = 1e-3

# Integrator params
if     space_time=="flat"   const initial_τs, final_τs = [-5σ, -5σ], [5σ, 5σ]
elseif space_time=="quench" const initial_τs, final_τs = [max(σ*1e-2, switching_function_center - 5σ),
                                                          max(σ*1e-2, switching_function_center - 5σ)], 
                                                         [switching_function_center + 5σ, 
                                                          switching_function_center + 5σ] end
const int_tol = 1e-3
const maxevals = 50000

# Complex contour params
const ε_contour = 1e-2
const deform_func_name = "cos4"