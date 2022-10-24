# WightmanFunction params
const space_time = "quench"
if     space_time=="quench" const with_derivative_coupling = true
elseif space_time=="flat"   const with_derivative_coupling = false end
# const ε_W = 1e-1

# Switching function params
const σ = 1.0
const switching_func_name = "gauss"

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
const ε_numeric_derivative = 1e-4

# Integrator params
if     space_time=="flat"   const initial_τs, final_τs = [-5σ, -5σ], [5σ, 5σ]
elseif space_time=="quench" const initial_τs, final_τs = [σ*1e-2, σ*1e-2], [5σ, 5σ] end
const int_tol = 1e-3

# Complex contour params
const ε_contour = 1e-1
const deform_func_name = "cos4"
