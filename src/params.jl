using DataFrames

calculating_concurrence = false
calculating_mutual_information = true

const experiment_name = "below_pulse"
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

# Detector frequencies 
nΩ = 20 # Number of frequencies
# Ω0, Ωf = calculating_mutual_information ? (-1/σ, 5/σ) : (10/σ, 20/σ)
Ω0, Ωf = (10/σ, 20/σ)
Ωs = LinRange(Ω0, Ωf, nΩ)

# Detector Initial positions
nχ = 20 #Number of initial conditions to iterate over
const χ0A = 2σ
const y0A = asinh(χ0A/b)
if χ0A == 0 println("χ0A is zero and that will create problems with rindler and quench"); throw(Exception) end
if     space_time in ["quench", "rindler"]  χ0B0, χ0Bf = abs(χ0A) + 1.05σ, abs(χ0A) + 1.5σ
elseif space_time == "flat_2d"              χ0B0, χ0Bf = abs(χ0A) + 1σ   , abs(χ0A) + 1.5σ end
χ0Bs =  LinRange(χ0B0, χ0Bf, nχ)
y0B0, y0Bf = asinh(χ0B0/b), asinh(χ0Bf/b)

# Switching functions parameters continued
const ηcA_or_τcA = occursin("above_pulse", experiment_name) ? 4y0A :
                   occursin("at_pulse"   , experiment_name) ? 1y0A :
                   occursin("below_pulse", experiment_name) ? 0.2865y0A : 1y0A
const ηcB_or_τcB = ηcA_or_τcA

# Coupling strength of detector
const λ = 1e-2

# Delta for numerical derivation 
const ε_numeric_derivative = 1e-3

# Integrator params
Δτ = switching_func_name == "cos4" ? 0.5σ : 5σ
const rtol = 1e-2
const maxevals = 500000

# Complex contour params
const ε_contour = 5e-4

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
                   "y0A"                  => [y0A],
                   "Ω0"                   => [Ω0], 
                   "Ωf"                   => [Ωf], 
                   "χ0B0"                 => [χ0B0], 
                   "χ0Bf"                 => [χ0Bf],
                   "y0B0"                 => [y0B0],
                   "y0Bf"                 => [y0Bf],
                   "λ"                    => [λ],
                   "Numeric_Derivative_ε" => [ε_numeric_derivative],
                #    "initial_τs"           => [initial_τs],
                #    "final_τs"             => [final_τs],
                   "rtol"                 => [rtol],
                   "maxevals"             => [maxevals],
                   "Contour_ε"            => [ε_contour])
