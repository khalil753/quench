using HCubature, Plots, QuantumInformation
include("tools/2_D&W.jl")
include("tools/3_LM_getters.jl")
include("tools/SwitchingFuncs.jl")
include("tools/DeformFuncs.jl")
include("params.jl")

integrate(f::Function) = hcubature(f, initial_τs, final_τs, maxevals=100000, rtol=int_tol)[1]

function initialize_stuff()
    Wττ′s, Dττ′ = create_distributions(_Ws[wightman_funct_name], _Ds[wightman_funct_name], 
                                       χ0A, χ0B, b, 
                                       numeric_derivative_ε)
    χ(τ) = χs[switching_func_name](τ/σ)
    df = deform_funcs[deform_func_name]
    dist_func = distance_funcs[dist_func_name]
    return Wττ′s, Dττ′, χ, df, dist_func
end

function main()
    Wττ′s, Dττ′, χ, df, dist_func = initialize_stuff()
    Cplify(l_or_m) = complexify_l_or_m(l_or_m, df, dist_func, ε_contour) # I'm writing this to shorten the name and to be able to use map_dict that requires a functions with a single argument

    Ωs = LinRange(1,5,10)
    Cs, ρs, = [], []
    for (i,Ω) in enumerate(Ωs)
        println("Doing Ω number $i: Ω = $Ω")
        m, ls = get_m(Dττ′, λ, Ω, χ), get_ls(Wττ′s, λ, Ω, χ)
        m, ls = Cplify(m)           , map_dict(Cplify, ls)
        M, Ls = integrate(m)        , map_dict(integrate, ls)

        ρ = [1 - Ls["AA"] - Ls["BB"]              0          0   conj(M);
                                   0       Ls["AA"]   Ls["AB"]         0;
                                   0  conj(Ls["AB"])  Ls["BB"]         0;
                                   M              0          0         0]

        my_concurrence(ρ) = 2*max(abs(ρ[1,4]) - real(ρ[2,2] + ρ[3,3])/2, 0)
        push!(Cs, my_concurrence(ρ))
        push!(ρs, ρ)
    end

    display(plot(Ωs, Cs))
    return ρs, Ωs, Cs
end

# if abspath(PROGRAM_FILE) == @__FILE__
ρ, Ωs, Cs = main()
# end
