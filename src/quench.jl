using QuadGK, HCubature, ProgressBars, Dates, CSV, CairoMakie, JLD
include("params.jl")
include("tools/3_LM_getters.jl")
include("tools/SwitchingFuncs.jl")
include("tools/Memoized_Integrator.jl")

println("\n", experiment_name, "\n")
if !save_plots println("Warning, you're running the simulation without saving the plots\n") end

integrate = MemoizedIntegrator(initial_τs, final_τs, maxevals, rtol)

χs = initialize_switching_funcs(switching_func_name, σ, switching_function_center_A, switching_function_center_B)  

ρs = []
Cs , Ns  = zeros(nΩ, nχ), zeros(nΩ, nχ)
PBs, PAs = zeros(nΩ, nχ), zeros(nΩ)
run_duration = begin
@elapsed for (i, Ω) in tqdm(enumerate(Ωs))
    for (j, χ0B) in tqdm(enumerate(χ0Bs))
        XA, XB = initialize_trajs(space_time, χ0A, χ0B, b)
        W = deform(_Ws[space_time], ε_contour)
        Ws, D = initialize_distributions(W, XA, XB, with_derivative_coupling, ε_numeric_derivative)

        m, ls = get_m(D, λ, Ω, χs), get_ls(Ws, λ, Ω, χs); delete!(ls, "AB")
        M, Ls = integrate(m), map_dict(integrate, ls)

        ρ = get_ρ(M, Ls)
        push!(ρs, ρ)
        Cs[i,j], Ns[i,j] = concurrence(ρ), negativity(ρ)
        PAs[i], PBs[i,j] = real(ρ[2,2]), real(ρ[3,3])
    end
end
end

path = "plots/new_plots/quench/"
img_name = make_img(χ0Bs, Ωs, Cs/λ^2, path, experiment_name, save_plots)
if save_data  save("plots/new_plots/quench/$(img_name).jld", "data", Cs) end
if save_plots store_in_df(path, "df.csv", params, [img_name], [run_duration]) end




