using DataFrames
const experiment_name = "new_paper_accelerated_detector_test"

# WightmanFunction params
const space_time = "flat"

const with_derivative_coupling = false 

# Switching function params
const switching_func_name = "gauss"
const σ = 1.0
const switching_function_center_A = 0.0σ
const switching_function_center_B = 0.0σ

# Detector frequencies and Initial positions
nΔL = 30 # Number of separations to iterate over
Ωs = [0.5/σ, 1.2/σ, 2.0/σ]
# ΔLss = Dict(0.01/σ => LinRange(0σ, 1.4σ, nΔL),
#              0.5/σ => LinRange(0σ,   2σ, nΔL),
#                2/σ => LinRange(0σ,   3σ, nΔL))
# C_ranges = Dict(0.01/σ => [0,  1.2],
#                  0.5/σ => [0,  0.6],
#                    2/σ => [0, 0.04])
ΔLss = Dict(0.5/σ => LinRange(0σ, 0.6σ, nΔL),
            1.2/σ => LinRange(0σ, 1.2σ, nΔL),
              2/σ => LinRange(0σ, 1.1σ, nΔL))
C_ranges = Dict(0.5/σ => [0,  3.0],
                1.2/σ => [0,  0.9],
                  2/σ => [0,  0.15])

# Ωs, ΔLs = [0.01/σ, 0.5/σ, 2.0/σ], LinRange(ΔL0, ΔLf, nχ)
nΩ = length(Ωs)
a = 3.0
χ0 = 1/a ; if χ0 == 0 println("χ0A is zero and that will create problems with rindler and quench"); throw(Exception) end

# Coupling strength of detector
const λ = 1e-1

# Delta for numerical derivation 
const ε_numeric_derivative = 5e-4

# Integrator params
Δτ = switching_func_name == "cos4" ? 0.5σ : 5σ
const initial_τs, final_τs = [switching_function_center_A - Δτ, switching_function_center_B - Δτ],
                             [switching_function_center_A + Δτ, switching_function_center_B + Δτ]
const rtol = 1e-3
const maxevals = 5000000

# Complex contour params
const ε_contour = 1e-3

params = DataFrame("Space_Time"           => [space_time],
                   "Derivative_Coupling"  => [with_derivative_coupling],
                   "Switching_Func"       => [switching_func_name],
                   "σ"                    => [σ],
                   "SF_Center_A"          => [switching_function_center_A],
                   "SF_Center_B"          => [switching_function_center_B],
                   #  "b"                    => [b],
                   #  "nΩ"                   => [20],
                   "nΔL"                  => [nΔL],
                   #  "χ0A"                  => [χ0A],
                   #  "Ω0"                   => [Ω0], 
                   #  "Ωf"                   => [Ωf], 
                   #  "χ0B0"                 => [χ0B0], 
                   #  "χ0Bf"                 => [χ0Bf],
                   "λ"                    => [λ],
                   "Numeric_Derivative_ε" => [ε_numeric_derivative],
                   "initial_τs"           => [initial_τs],
                   "final_τs"             => [final_τs],
                   "rtol"                 => [rtol],
                   "maxevals"             => [maxevals],
                   "Contour_ε"            => [ε_contour])
