using HCubature, ProgressBars, Dates, CSV, CairoMakie, SpecialFunctions
include("../rindler_main/params2.jl")
include("../tools/3_LM_getters.jl")
include("../tools/SwitchingFuncs.jl")
include("../tools/Memoized_Integrator.jl")

integrate = MemoizedIntegrator(initial_τs, final_τs, maxevals, rtol)
χs = initialize_switching_funcs(switching_func_name, σ, switching_function_center_A, switching_function_center_B)  

ρs, run_durations, Cs = [], [], Dict()
for (i, Ω) in tqdm(enumerate(Ωs))
    ΔLs = ΔLss[Ω]
    Cs[Ω] = []
    run_duration = begin
    @elapsed  for (j, ΔL) in tqdm(enumerate(ΔLs))
        XA, XB = AcceleratedTrajectory(χ0, 0, 0, 0), AcceleratedTrajectory(χ0, ΔL, 0, 0)
        Ws, D  = initialize_distributions(_Ws[space_time], _Ds[space_time], XA, XB)
        if with_derivative_coupling Ws, D = add_crossed_derivatives(Ws, D, ε_numeric_derivative) end

        m, ls = get_m(D, λ, Ω, χs, ε_contour), get_ls(Ws, λ, Ω, χs, ε_contour)
        delete!(ls, "AB")
        M, Ls = integrate(m)                 , map_dict(integrate, ls)

        ρ = get_ρ(M, Ls)
        push!(ρs, ρ)
        push!(Cs[Ω], concurrence(ρ))
    end
    end
    push!(run_durations, run_duration)
end

path = "plots/new_plots/rindler_plots"
img_names = plot_C_vs_L(path, experiment_name, ΔLss, Ωs, Cs)
store_in_df(path, "df.csv", params, img_names, run_durations)