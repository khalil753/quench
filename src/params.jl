using DataFrames

const experiment_name = "at_quench"
const save_plots = true
const save_data  = true

# WightmanFunction params
const space_time = "quench"

if     space_time == "quench"  const with_derivative_coupling = true
elseif space_time == "rindler" const with_derivative_coupling = true 
elseif space_time == "flat_2d" const with_derivative_coupling = false end

# Switching function params
const switching_func_name = "cos4"
const σ = 1.0
const using_ηs = true

# Quench regulator
const b = space_time=="quench" ? (5e-1)*σ : 0.0

# Detector frequencies and Initial positions
nΩ = 20 # Number of frquencies/initial conditions to iterate over
nχ = 20
const χ0A = 2σ 
const y0A = asinh(χ0A/b)
if χ0A == 0 println("χ0A is zero and that will create problems with rindler and quench"); throw(Exception) end
if     space_time in ["quench", "rindler"]  Ω0, Ωf, χ0B0, χ0Bf = 10/σ, 20/σ, abs(χ0A) + 1.05σ, abs(χ0A) + 1.5σ
elseif space_time == "flat_2d"              Ω0, Ωf, χ0B0, χ0Bf = 10/σ, 20/σ, abs(χ0A) + 1σ   , abs(χ0A) + 1.5σ end
Ωs, χ0Bs = LinRange(Ω0, Ωf, nΩ), LinRange(χ0B0, χ0Bf, nχ)

# Switching functions parameters continued
const ηcA_or_τcA = 1.0y0A
const ηcB_or_τcB = 1.0y0A

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
