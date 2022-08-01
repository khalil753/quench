# WightmanFunction params
wightman_funct_name = "quench"
dist_ε = 1e-1

# Switching function params
σ = 0.5
switching_func_name = "cos"

# Quench regulator
b = 1.0

# Initial detector positions
χ0A, χ0B = 1.0, 5.0

# Detector frequencies
Ω = 1.5
# ΩA, ΩB = 1.0, 1.0

# Coupling strength of detector
λ = 0.1

# Delta for numerical derivation 
numeric_derivative_ε = 1e-5

# Integrator params
initial_τs, final_τs =  [0.1,0.1], [σ, σ]
int_tol = 1e-4
