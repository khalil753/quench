# WightmanFunction params
const space_time = "quench"
if     space_time=="quench" const with_derivative_coupling = true
elseif space_time=="flat"   const with_derivative_coupling = false end
# const ε_W = 1e-1

# Switching function params
const σ = 1.0
const switching_func_name = "cos4"
const switching_function_center = 2σ
# const switching_function_center_A = 0σ
# const switching_function_center_B = 0σ

# Quench regulator
const b = 0.01

# Initial detector positions
const χ0A, χ0B = 1.0, 2.0

# Detector frequencies
const Ω = 1.5
# ΩA, ΩB = 1.0, 1.0

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

# Complex contour params
const ε_contour = 1e-1
const deform_func_name = "cos4"