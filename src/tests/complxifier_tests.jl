using HCubature, Plots, QuantumInformation
include("../params.jl")
include("../tools/misc.jl")
include("../tools/1_Trajectories.jl")
include("../tools/2_D&W.jl")
include("../tools/3_LM_getters.jl")
include("../tools/SwitchingFuncs.jl")
include("../tools/DeformFuncs.jl")

f(x) = _gauss(x[1])*_gauss(x[2])/(x[1] - x[2])^2
pole_distance(x) = (x[1] - x[2])^2
f_c = complexify(f, deform_funcs["cos4"], pole_distance, 1e-2)

n = 5
times = []
for τ1 in 0:0.1:1
    for τ2 in 0:0.1:1
        push!(times, @btime f_c([$τ1, $τ2]))
    end
end

