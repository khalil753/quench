using HCubature, Plots, QuantumInformation
include("tools/misc.jl")
include("tools/params.jl")
include("tools/D&W.jl")
include("tools/Trajectories.jl")
include("tools/SwitchingFuncs.jl")

function create_Ws(χ0A, χ0B, b)
    Ws = Dict()
    Ws["AB"] = WightmanFunction(χ0A, χ0B, b)
    Ws["AA"] = WightmanFunction(χ0A, χ0A, b)
    Ws["BB"] = WightmanFunction(χ0B, χ0B, b)
    return Ws
end

integrate(f::Function) = hcubature(f, initial_τs, final_τs, rtol=int_tol)[1]
get_crossed_derivative(f) = get_crossed_derivative(f, ε)

χ = χs[switching_func_name]
Ωs = LinRange(1,5,10)
Cs, ρs, Lss, Ms = [], [], [], []
for Ω in Ωs
    Ws = create_Ws(χ0A, χ0B, b)
    Wττ′s = map_dict(get_crossed_derivative, Ws)
    Dττ′  = get_crossed_derivative(Propagator(χ0A, χ0B, b))

    ls = get_ls(Wττ′s, λ, Ω, χ)
    m  = get_m(Dττ′, λ, Ω, χ)

    M = integrate(m)
    Ls = map_dict(integrate, ls)

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