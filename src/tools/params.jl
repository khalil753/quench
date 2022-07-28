# Switching function params
σ = 0.5
switching_func_name = "cos"

# Quench regulator
b = 1.0

# Initial detector positions
χ0A, χ0B = 1.0, 5.0

# Detector frequencies
Ω = 1.0
# ΩA, ΩB = 1.0, 1.0

# Coupling strength of detector
λ = 0.1

# Delta for numerical derivation 
ε = 1e-5

# Integrator params
initial_τs, final_τs =  [0.1,0.1], [σ, σ]
int_tol = 1e-4
