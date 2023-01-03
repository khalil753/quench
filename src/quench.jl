using HCubature, Plots, ProgressBars, DataFrames#, QuantumInformation
include("params.jl")
# include("tools/misc.jl")
# include("tools/1_Trajectories.jl")
# include("tools/2_D&W.jl")
include("tools/3_LM_getters.jl")
include("tools/SwitchingFuncs.jl")

function get_trajs(space_time, χ0A, χ0B, b)
    if     space_time=="quench"  XA, XB = QuenchTrajectory(χ0A, b)   , QuenchTrajectory(χ0B, b)
    elseif space_time=="flat"    XA, XB = InertialTrajectory(0, 0, 0), InertialTrajectory(χ0B, 0, 0) 
    elseif space_time=="rindler" XA, XB = AcceleratedTrajectory(χ0A,0,0), AcceleratedTrajectory(χ0B,0,0) end
    return XA, XB
end

χs = initialize_switching_funcs(switching_func_name, σ, switching_function_center_A, switching_function_center_B)  
integrate(f) = hcubature(f, initial_τs, final_τs, maxevals=maxevals, rtol=int_tol)[1]

ρs = []
Cs, Ns, = zeros(length(Ωs), length(χ0Bs)), zeros(length(Ωs), length(χ0Bs))
PAs, PBs = zeros(length(Ωs)), zeros(length(Ωs), length(χ0Bs))
for (i, Ω) in tqdm(enumerate(Ωs))
    for (j, χ0B) in tqdm(enumerate(χ0Bs))
        XA, XB = get_trajs(space_time, χ0A, χ0B, b)
        Ws, D  = initialize_distributions(_Ws[space_time], _Ds[space_time], XA, XB)
        if with_derivative_coupling Ws, D = add_crossed_derivatives(Ws, D, ε_numeric_derivative) end

        m, ls = get_m(D, λ, Ω, χs, ε_contour), get_ls(Ws, λ, Ω, χs, ε_contour)
        M, Ls = integrate(m)                 , map_dict(integrate, ls)

#                                 00        01              10        11
        ρ = [1 - Ls["AA"] - Ls["BB"]         0               0   conj(M); #00
                                   0  Ls["BB"]   conj(Ls["AB"])        0; #01
                                   0  Ls["AB"]        Ls["AA"]         0; #10
                                   M         0               0         0] #11

        Cs[i,j] = concurrence(ρ)
        Ns[i,j] = negativity(ρ)
        push!(ρs, ρ)
        PAs[i]   = real(ρ[2,2])
        PBs[i,j] = real(ρ[3,3])
    end
end

p = contourf(χ0Bs, Ωs, Cs, linewidth=-0.0, xlabel="χB₀", ylabel="Ω", title="Concurrence")
display(p)
savefig(p, "plots/$(space_time)_concurrence_vs_Ω_χ.png")


