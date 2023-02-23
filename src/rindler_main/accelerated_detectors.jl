using HCubature

using HCubature, ProgressBars, Dates, CSV, CairoMakie, SpecialFunctions
include("../rindler_main/params2.jl")
include("../tools/3_LM_getters.jl")
include("../tools/SwitchingFuncs.jl")
include("../tools/Memoized_Integrator.jl")

integrate = MemoizedIntegrator(initial_τs, final_τs, maxevals, rtol)
χs = initialize_switching_funcs(switching_func_name, σ, switching_function_center_A, switching_function_center_B)  

ρs, run_durations, Cs = [], [], Dict()

for (i, Ω) in tqdm(enumerate(Ωs))
    function p1(ss)
        s = ss[1]
        β = 2Ω/a
        α = 1/(a*σ^2)
        λ^2*a*σ/4π^(3/2)*(cos(s*β)*exp(-s^2*α)*(sinh(s)^2 - s^2)/(s^2*sinh(s)^2))
    end
    P = hcubature(p1, [0.001], [5σ], maxevals=maxevals, rtol=rtol)[1] + 
    λ^2/4π*(exp(-Ω^2*σ^2) - √π*Ω*σ*erfc(Ω*σ))

    ΔLs = ΔLss[Ω]
    Cs[Ω] = []
    run_duration = begin
    @elapsed for (j, ΔL) in tqdm(enumerate(ΔLs))
        function D_parallel(x,y) 
            a/32π^2(((a*ΔL/2 - exp(-x*a/2)*sinh(y*a/2))*(a*ΔL/2 + exp(x*a/2)*sinh(y*a/2)) - im*ε_contour)^(-1) + 
                    ((a*ΔL/2 + exp(-x*a/2)*sinh(y*a/2))*(a*ΔL/2 - exp(x*a/2)*sinh(y*a/2)) - im*ε_contour)^(-1)) 
        end        
        function F(xs)
            x,y = xs
            exp(-(x^2 + y^2)/4σ^2 - im*x*Ω)*D_parallel(x,y)    
        end
        ρ14 = -λ^2 * hcubature(F, [-5σ, 1], [5σ, 5σ], maxevals=maxevals, rtol=rtol)[1]
        println("ρ14 = $ρ14")
        ρ = get_ρ(ρ14, Dict("AA" => P, "BB" => P))
        push!(ρs, ρ)
        push!(Cs[Ω], concurrence(ρ))
    end
    end
    push!(run_durations, run_duration)
end

path = "new_plots/rindler_plots"
img_names = plot_C_vs_L(path, experiment_name, ΔLss, Ωs, Cs)
store_in_df(path, "df.csv", params, img_names, run_durations)



