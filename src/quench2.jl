using HCubature, Plots, ProgressBars, Dates, CSV
include("params2.jl")
include("tools/3_LM_getters.jl")
include("tools/SwitchingFuncs.jl")
include("tools/Memoized_Integrator.jl")

integrate = MemoizedIntegrator(initial_τs, final_τs, maxevals, int_tol)
χs = initialize_switching_funcs(switching_func_name, σ, switching_function_center_A, switching_function_center_B)  

ρs = []
Cs, Ns, = zeros(nΩ, nχ), zeros(nΩ, nχ)
PAs, PBs = zeros(length(Ωs)), zeros(length(Ωs), length(χ0Bs))
run_duration = begin
@elapsed for (i, Ω) in tqdm(enumerate(Ωs))
    for (j, χ0B) in tqdm(enumerate(χ0Bs))
        XA, XB = initialize_trajs(space_time, χ0A, χ0B, b)
        Ws, D  = initialize_distributions(_Ws[space_time], _Ds[space_time], XA, XB)
        if with_derivative_coupling Ws, D = add_crossed_derivatives(Ws, D, ε_numeric_derivative) end

        m, ls = get_m(D, λ, Ω, χs, ε_contour), get_ls(Ws, λ, Ω, χs, ε_contour) 
        M, Ls = integrate(m)                 , map_dict(integrate, ls)

        ρ = get_ρ(M, Ls)
        push!(ρs, ρ)
        Cs[i,j], Ns[i,j] = concurrence(ρ), negativity(ρ)
        PAs[i], PBs[i,j] = real(ρ[2,2]), real(ρ[3,3])
    end
end
end

display(Cs)
img_name = make_img(χ0Bs, Ωs, Cs)
store_in_df(params, img_name, run_duration)




