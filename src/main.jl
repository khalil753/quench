using QuadGK, HCubature, ProgressBars, Dates, CSV, CairoMakie, JLD, Printf
include("params.jl")
include("tools/3_LM_getters.jl")
include("tools/SwitchingFuncs.jl")
include("tools/Memoized_Integrator.jl")

println("\n", experiment_name, "\n")
if !save_plots println("Warning, you're running the simulation without saving the plots\n") end

integrate = MemoizedIntegrator(maxevals, rtol)

ρs = []
Cs , Ns, MIs  = zeros(nΩ, nχ), zeros(nΩ, nχ), zeros(nΩ, nχ)
PBs, PAs = zeros(nΩ, nχ), zeros(nΩ)
run_duration = begin
@elapsed for (i, Ω) in tqdm(enumerate(Ωs))
    for (j, χ0B) in tqdm(enumerate(χ0Bs))
        initial_τs, final_τs = initialize_integration_ranges(ηcA_or_τcA, ηcB_or_τcB, χ0A, χ0B, b, using_ηs, space_time, switching_func_name)
        χs = initialize_switching_funcs(switching_func_name, σ, ηcA_or_τcA, ηcB_or_τcB, using_ηs, χ0A, χ0B)

        XA, XB = initialize_trajs(space_time, χ0A, χ0B, b)
        W = deform(_Ws[space_time], ε_contour)
        Ws, D = initialize_distributions(W, XA, XB, with_derivative_coupling, ε_numeric_derivative)

        m, ls = get_m(D, λ, Ω, χs), get_ls(Ws, λ, Ω, χs); delete!(ls, "AB")
        M, Ls = integrate_m_and_ls(m, ls, initial_τs, final_τs, integrate) 

        ρ = get_ρ(M, Ls)
        push!(ρs, ρ)
        Cs[i,j], Ns[i,j], MIs[i,j] = concurrence(ρ), negativity(ρ), mutual_information(ρ)
        PAs[i], PBs[i,j] = real(ρ[2,2]), real(ρ[3,3])
    end
end
end

plot_folder = "quench_plots"
path = "plots/$plot_folder"

if calculating_concurrence  
    println("Making Concurrence Plot...")      
    img_name1 = make_img(χ0Bs, Ωs, Cs/λ^2 + 0.0000000001*rand(size(Cs)...), path, "concurrence_"*experiment_name       , save_plots)
    if save_data save("plots/$plot_folder/$(img_name1).jld", "data", Cs ) end
    if save_plots store_in_df(path, "df.csv", params, [img_name1], [run_duration]) end
end
if calculating_mutual_information 
    println("Making Mutual Information Plot...")      
    img_name2 = make_img(χ0Bs, Ωs, MIs, path, "mutual_information_"*experiment_name, save_plots)
    if save_data save("plots/$plot_folder/$(img_name2).jld", "data", MIs) end
    if save_plots store_in_df(path, "df.csv", params, [img_name2], [run_duration]) end
end




