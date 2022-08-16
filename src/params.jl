# WightmanFunction params
const space_time = "quench"
# const ε_W = 1e-1

# Switching function params
const σ = 1.0
const switching_func_name = "cos4"

# Quench regulator
const b = 0.1

# Initial detector positions
const χ0A, χ0B = 1.0, 15.0

# Detector frequencies
const Ω = 1.5
# ΩA, ΩB = 1.0, 1.0

# Coupling strength of detector
const λ = 1e-1

# Delta for numerical derivation 
const ε_numeric_derivative = 1e-4

# Integrator params
const initial_τs, final_τs =  [0.1,0.1], [σ, σ]
const int_tol = 1e-3

# Complex contour params
const ε_contour = 1e-2
const deform_func_name = "cos4"
