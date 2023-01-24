using HCubature, ProgressBars, Dates, CSV, CairoMakie, SpecialFunctions
include("../rindler_main/params2.jl")
include("../tools/3_LM_getters.jl")
include("../tools/SwitchingFuncs.jl")
include("../tools/Memoized_Integrator.jl")

integrate = MemoizedIntegrator(initial_τs, final_τs, maxevals, rtol)
χs = initialize_switching_funcs(switching_func_name, σ, switching_function_center_A, switching_function_center_B)  

ρs = []
Cs  , Ns   = zeros(nΩ, nχ), zeros(nΩ, nχ)
Cs_2, Ns_2 = zeros(nΩ, nχ), zeros(nΩ, nχ)
PBs , PAs  = zeros(nΩ, nχ), zeros(nΩ)
run_duration = begin
@elapsed for (i, Ω) in tqdm(enumerate(Ωs))
    ΔLs = ΔLss[Ω]
    # a = 1/χ0
    # function p1(ss)
    #     s = ss[1]
    #     β = 2Ω/a
    #     α = 1/(a*σ^2)
    #     λ^2*a*σ/4π^(3/2)*(cos(s*β)*exp(-s^2*α)*(sinh(s)^2 - s^2)/(s^2*sinh(s)^2))
    # end
    # P = hcubature(p1, [0.001], [5σ], maxevals=10000000, rtol=1e-10)[1] + 
    #     λ^2/4π*(exp(-Ω^2*σ^2) - √π*Ω*σ*erfc(Ω*σ))
    # println("P = $P")
    for (j, ΔL) in tqdm(enumerate(ΔLs))
        XA, XB = AcceleratedTrajectory(χ0, 0, 0, 0), AcceleratedTrajectory(χ0, ΔL, 0, 0)
        Ws, D  = initialize_distributions(_Ws[space_time], _Ds[space_time], XA, XB)
        if with_derivative_coupling Ws, D = add_crossed_derivatives(Ws, D, ε_numeric_derivative) end

        m, ls = get_m(D, λ, Ω, χs, ε_contour), get_ls(Ws, λ, Ω, χs, ε_contour)
        M, Ls = integrate(m)                 , map_dict(integrate, ls)

        ρ = get_ρ(M, Ls)
        push!(ρs, ρ)
        Cs[i,j], Ns[i,j] = concurrence(ρ), negativity(ρ)
        PAs[i], PBs[i,j] = real(ρ[2,2]), real(ρ[3,3])
        
        # function D_parallel(x,y) 
        #     a/32π^2(((a*ΔL/2 - exp(-x*a/2)*sinh(y*a/2))*(a*ΔL/2 + exp(x*a/2)*sinh(y*a/2)) - im*1e-1)^(-1) + 
        #             ((a*ΔL/2 + exp(-x*a/2)*sinh(y*a/2))*(a*ΔL/2 - exp(x*a/2)*sinh(y*a/2)) - im*1e-1)^(-1)) 
        # end        
        # function F(xs)
        #     x,y = xs
        #     exp(-(x^2 + y^2)/4σ^2 - im*x*Ω)*D_parallel(x,y)    
        # end
        # ρ14 = -λ^2 * hcubature(F, [0.01, -5σ], [5σ, 5σ], maxevals=100000000, rtol=1e-10)[1]
        # println("ρ14 = $ρ14")
        # ρ = get_ρ(ρ14, Dict("AA" => P, "BB" => P))
        # Cs_2[i,j], Ns_2[i,j] = concurrence(ρ), negativity(ρ)
    end
end
end

path = "new_plots/rindler_plots"
img_name   = plot_C_vs_L(path, ΔLs, Ωs, Cs')
# img_name_2 = plot_C_vs_L(ΔLs, Ωs, Cs_2')
store_in_df(path, "df.csv", params, img_name, run_duration)
# store_in_df(path, "df.csv", params, img_name_2, run_duration)



