using HCubature, Plots, QuantumInformation
include("tools/misc.jl")
include("tools/params.jl")
include("tools/D&W.jl")
include("tools/Trajectories.jl")
include("tools/SwitchingFuncs.jl")

X = InertialTrajectory(0.0, 0.0, 0.0)
W = WightmanFunction(X, X)

l = get_l(W, 1, 1, χs["one"])

eps = 1e-6
f(t) = exp(-im*t)/(4π*(t - im*eps)^2)
f_vec(t) = f(t[1])

t = LinRange(0,1e-5,1000)
g(t) = l([t,0.0])
plot(t, [abs ∘ g,abs ∘ f])