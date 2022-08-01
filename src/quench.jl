using HCubature, Plots, QuantumInformation
include("tools/misc.jl")
include("tools/params.jl")
include("tools/D&W.jl")
include("tools/Trajectories.jl")
include("tools/SwitchingFuncs.jl")

integrate(f::Function) = hcubature(f, initial_τs, final_τs, rtol=int_tol)[1]

Wττ′s, Dττ′ = create_distributions(_Ws, _Ds, wightman_funct_name, 
                                   χ0A, χ0B, b, 
                                   dist_ε, numeric_derivative_ε)

χ = χs[switching_func_name]
Ωs = LinRange(1,5,10)
Cs, ρs, Lss, Ms = [], [], [], []
for Ω in Ωs
    m, ls = get_m(Dττ′, λ, Ω, χ), get_ls(Wττ′s, λ, Ω, χ)
    M, Ls = integrate(m)        , map_dict(integrate, ls)

    ρ = [1 - Ls["AA"] - Ls["BB"]              0          0   conj(M);
                               0       Ls["AA"]   Ls["AB"]         0;
                               0  conj(Ls["AB"])  Ls["BB"]         0;
                               M              0          0         0]
    for i in size(ρ)[1]; ρ[i,i] = real(ρ[i,i,]) end

    push!(Cs, concurrence(ρ))
    push!(ρs, ρ)
    push!(Lss, Ls)
    push!(Ms, M)
end

plot(Ωs, Cs)