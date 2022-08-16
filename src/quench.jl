using HCubature, Plots#, QuantumInformation
include("params.jl")
include("tools/misc.jl")
include("tools/1_Trajectories.jl")
include("tools/2_D&W.jl")
include("tools/3_LM_getters.jl")
include("tools/SwitchingFuncs.jl")
include("tools/DeformFuncs.jl")

Wττ′s, Dττ′, χ, df, dist_func = initialize_stuff()

integrate(f) = hcubature(f, initial_τs, final_τs, maxevals=100000, rtol=int_tol)[1]
# Cplify(l_or_m) = complexify_l_or_m(l_or_m, df, dist_func, ε_contour) # I'm writing this to shorten the name and to be able to use map_dict that requires a functions with a single argument
Cplify(l_or_m) = complexify_l_or_m(l_or_m, ε_contour) # I'm writing this to shorten the name and to be able to use map_dict that requires a functions with a single argument

Ωs = LinRange(0,1,7)
Ns, Cs, ρs, lss, ms = [], [], [], [], []
for (i,Ω) in enumerate(Ωs)
    println("Doing Ω number $i: Ω = $Ω")
    m, ls = get_m(Dττ′, λ, Ω, χ), get_ls(Wττ′s, λ, Ω, χ)
    m, ls = Cplify(m)           , map_dict(Cplify, ls)
    push!(ms,m); push!(lss,ls)
    M, Ls = integrate(m)        , map_dict(integrate, ls)

    ρ = [1 - Ls["AA"] - Ls["BB"]              0          0   conj(M);
                               0       Ls["AA"]   Ls["AB"]         0;
                               0  conj(Ls["AB"])  Ls["BB"]         0;
                               M              0          0         0]

    push!(Cs, concurrence(ρ))
    push!(Ns, negativity(ρ))
    push!(ρs, ρ)
end

p = plot(Ωs, Cs)
display(p)
savefig(p, "plots/concurrence_vs_Ω.png")


