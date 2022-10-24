using HCubature, Plots, ProgressBars#, QuantumInformation
include("params.jl")
# include("tools/misc.jl")
# include("tools/1_Trajectories.jl")
# include("tools/2_D&W.jl")
include("tools/3_LM_getters.jl")
include("tools/SwitchingFuncs.jl")
include("tools/DeformFuncs.jl")

χ(τ) = χs[switching_func_name](τ/σ)
df = deform_funcs[deform_func_name]
dist_func = distance_funcs[space_time]

integrate(f) = hcubature(f, initial_τs, final_τs, maxevals=20000, rtol=int_tol)[1]
# Cplify(l_or_m) = complexify_l_or_m(l_or_m, df, dist_func, ε_contour) # I'm writing this to shorten the name and to be able to use map_dict that requires a functions with a single argument
Cplify(l_or_m) = complexify_l_or_m(l_or_m, ε_contour) # I'm writing this to shorten the name and to be able to use map_dict that requires a functions with a single argument

ρs = []
if     space_time=="quench" Ωs, χ0Bs = LinRange(-1/σ, 10/σ, 20), LinRange(χ0A + 0.5σ, χ0A + 5σ, 20)
elseif space_time=="flat"   Ωs, χ0Bs = LinRange(-3/σ,  3/σ, 20), LinRange(0.5σ, 2σ, 20) end
Cs, Ns = zeros(length(Ωs), length(χ0Bs)), zeros(length(Ωs), length(χ0Bs))
for (i, Ω) in tqdm(enumerate(Ωs))
    for (j, χ0B) in tqdm(enumerate(χ0Bs))
        if     space_time=="quench" XA, XB = QuenchTrajectory(χ0A, b)   , QuenchTrajectory(χ0B, b)
        elseif space_time=="flat"   XA, XB = InertialTrajectory(0, 0, 0), InertialTrajectory(χ0B, 0, 0) end
        Ws, D = initialize_distributions(_Ws[space_time], _Ds[space_time], XA, XB)
        if with_derivative_coupling Ws, D = add_crossed_derivatives(Ws, D, ε_numeric_derivative) end

        m, ls = get_m(D, λ, Ω, χ)   , get_ls(Ws, λ, Ω, χ)
        m, ls = Cplify(m)           , map_dict(Cplify, ls)
        M, Ls = integrate(m)        , map_dict(integrate, ls)

        ρ = [1 - Ls["AA"] - Ls["BB"]              0          0   conj(M);
                                   0       Ls["AA"]   Ls["AB"]         0;
                                   0  conj(Ls["AB"])  Ls["BB"]         0;
                                   M              0          0         0]

        Cs[i,j] = concurrence(ρ)
        Ns[i,j] = negativity(ρ)
        push!(ρs, ρ)
    end
end

p = contourf(Ωs, χ0Bs, Cs', linewidth=-0.0)
display(p)
savefig(p, "plots/$(space_time)_concurrence_vs_Ω_χ.png")


