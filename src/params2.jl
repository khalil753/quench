using DataFrames

# WightmanFunction params
const space_time = "rindler"

with_derivative_coupling = false 

# Switching function params
const switching_func_name = "gauss"
const σ = 1.0
const switching_function_center_A = 0.5σ
const switching_function_center_B = 0.5σ

# Detector frequencies and Initial positions
nΩ, nχ = 20, 20 # Number of frquencies/initial conditions to iterate over
const χ0A = 50σ 
if χ0A == 0 println("χ0A is zero and that will create problems with rindler and quench"); throw(Exception) end
Ω0, Ωf = -1/σ, 30/σ
χ0B0, χ0Bf = χ0A + 0.5σ, χ0A + 1.5σ
Ωs, χ0Bs = LinRange(Ω0, Ωf, nΩ), LinRange(χ0B0, χ0Bf, nχ)

# Coupling strength of detector
const λ = 1e-2

# Delta for numerical derivation 
const ε_numeric_derivative = 1e-3

# Integrator params
Δτ = switching_func_name == "cos4" ? 0.5σ : 5σ
const initial_τs, final_τs = [switching_function_center_A - Δτ, switching_function_center_B - Δτ],
                                                      [switching_function_center_A + Δτ, switching_function_center_B + Δτ]
const int_tol = 1e-4
const maxevals = 500000

# Complex contour params
const ε_contour = 1e-2

params = DataFrame("Space_Time"           => [space_time],
                   "Derivative_Coupling"  => [with_derivative_coupling],
                   "Switching_Func"       => [switching_func_name],
                   "σ"                    => [σ],
                   "SF_Center_A"          => [switching_function_center_A],
                   "SF_Center_B"          => [switching_function_center_B],
                #    "b"                    => [b],
                   "nΩ"                   => [20],
                   "nχ"                   => [20],
                   "χ0A"                  => [χ0A],
                   "Ω0"                   => [Ω0], 
                   "Ωf"                   => [Ωf], 
                   "χ0B0"                 => [χ0B0], 
                   "χ0Bf"                 => [χ0Bf],
                   "λ"                    => [λ],
                   "Numeric_Derivative_ε" => [ε_numeric_derivative],
                   "initial_τs"           => [initial_τs],
                   "final_τs"             => [final_τs],
                   "int_tol"              => [int_tol],
                   "maxevals"             => [maxevals],
                   "Contour_ε"            => [ε_contour])
