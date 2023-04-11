using DataFrames

const experiment_name = "quench_near_pulse_activation_A_center=3σ"
const save_plots = false
const save_data  = false

# WightmanFunction params
const space_time = "quench"

if     space_time == "quench"  const with_derivative_coupling = true
elseif space_time == "rindler" const with_derivative_coupling = true 
elseif space_time == "flat"    const with_derivative_coupling = false end

# Switching function params
const switching_func_name = "cos4"
const σ = 1.0
const ηcA_or_τcA = 1.0σ
const ηcB_or_τcB = 1.0σ
const using_ηs = false

# Quench regulator
const b = (1e-1)*σ

# Detector frequencies and Initial positions
nΩ = 10 # Number of frquencies/initial conditions to iterate over
nχ = 12
const χ0A = 1σ 
if χ0A == 0 println("χ0A is zero and that will create problems with rindler and quench"); throw(Exception) end
if     space_time in ["quench", "rindler"]  Ω0, Ωf, χ0B0, χ0Bf = -1/σ, 30/σ, χ0A + 0.5σ, χ0A + 1.5σ
elseif space_time == "flat"                 Ω0, Ωf, χ0B0, χ0Bf = -1/σ, 30/σ, χ0A + 0.5σ, χ0A + 1.5σ end
Ωs, χ0Bs = LinRange(Ω0, Ωf, nΩ), LinRange(χ0B0, χ0Bf, nχ)

# Coupling strength of detector
const λ = 1e-2

# Delta for numerical derivation 
const ε_numeric_derivative = 1e-3

# Integrator params
Δτ = switching_func_name == "cos4" ? 0.5σ : 5σ
const rtol = 1e-2
const maxevals = 500000

# Complex contour params
const ε_contour = 5e-3

params = DataFrame("Space_Time"           => [space_time],
                   "Derivative_Coupling"  => [with_derivative_coupling],
                   "Switching_Func"       => [switching_func_name],
                   "σ"                    => [σ],
                   "SF_Center_A"          => [ηcA_or_τcA],
                   "SF_Center_B"          => [ηcB_or_τcB],
                   "b"                    => [b],
                   "nΩ"                   => [nΩ],
                   "nχ"                   => [nχ],
                   "χ0A"                  => [χ0A],
                   "Ω0"                   => [Ω0], 
                   "Ωf"                   => [Ωf], 
                   "χ0B0"                 => [χ0B0], 
                   "χ0Bf"                 => [χ0Bf],
                   "λ"                    => [λ],
                   "Numeric_Derivative_ε" => [ε_numeric_derivative],
                #    "initial_τs"           => [initial_τs],
                #    "final_τs"             => [final_τs],
                   "rtol"                 => [rtol],
                   "maxevals"             => [maxevals],
                   "Contour_ε"            => [ε_contour])
