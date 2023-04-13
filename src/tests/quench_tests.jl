using HCubature, CairoMakie, SpecialFunctions, ProgressBars, Plots
include("../params.jl")
include("../tools/3_LM_getters.jl")
include("../tools/SwitchingFuncs.jl")
include("../tools/DeformFuncs.jl")
include("../tools/Memoized_Integrator.jl")

function get_P_Minkowski()
    """
    In this test I calculate the transition probability as a function of Ω, of an inertial observer
    in flat spacetime with a gaussian switching function.
    """
    λ = 1.0
    d = 5*σ    
    df = deform_funcs["cos2"]
    ε_contour = 1e-2

    initial_τs, final_τs =  [-d, -d], [d, d]
    integrate(f::Function) = hcubature(f, initial_τs, final_τs, maxevals=100000, rtol=rtol)[1]

    XA, XB = initialize_trajs("flat", 0, 0, b)
    W = deform(_Ws["flat"], ε_contour)
    Ws, _ = initialize_distributions(W, XA, XB, false, ε_numeric_derivative)    
    W = Ws["AA"]
    
    χ(τ) = switching_funcs["gauss"](τ/σ)

    Ωs = LinRange(0, 2.5, 30)
    P_Ms = []
    for (i, Ω) in enumerate(Ωs)
        println("\rDoing Ω number $i")
        l = get_l(W, λ, Ω, χ)
        P_M = integrate(l)
        push!(P_Ms, P_M)
    end
    x = max(imag.(P_Ms)...)
    println("$x")

    P(Ω) = λ^2/4π*(exp(-σ^2*Ω^2) - √π*σ*Ω*erfc(σ*Ω))

    p = Plots.plot(Ωs, [P.(Ωs), real.(P_Ms), imag.(P_Ms)], labels=["Theoretical" "Real part of numerical" "Imaginary part of numerical"])
    display(p)
    savefig(p, raw"plots\\test_plots\\\P_vs_Ω_Minkowski.png")
end

function get_M_vs_Ω_Minkowski()
    """
    In this test I calculate the transition probability as a function of Ω, of an inertial observer
    in flat spacetime with a gaussian switching function.
    """
    σ = 1.0
    d = 5*σ    
    df = deform_funcs["cos2"]
    ε_contour = 1e-2
    L = 10.0

    initial_τs, final_τs =  [-d, -d], [d, d]
    integrate(f) = hcubature(f, initial_τs, final_τs, maxevals=100000 , rtol=rtol)[1]

    XA, XB = initialize_trajs("flat", 0, L, b)
    W = deform(_Ws["flat"], ε_contour)
    _, D = initialize_distributions(W, XA, XB, false, ε_numeric_derivative)    
    χ(τ) = switching_funcs["gauss"](τ/σ)

    Ωs = LinRange(0, 2.5, 30)
    Ms = []
    for (i, Ω) in enumerate(Ωs)
        println("\rDoing Ω number $i: Ω = $Ω")
        m = get_m(D, λ, Ω, χ)
        M = integrate(m)
        push!(Ms, M)
    end

    M(Ω) = im*(λ^2)*σ/(4*√π*L)*exp(-(σ*Ω)^2 - L^2/(4*σ^2))*(erf(im*L/(2σ)) - 1)

    p = Plots.plot(Ωs, [(abs ∘ M).(Ωs), abs.(Ms)], labels=["theoretical" "numerical"])
    display(p)
    savefig(p, raw"plots\\test_plots\\M_vs_Ω_Minkowski.png")
    return Ωs, Ms, M.(Ωs)
end

function get_M_vs_L_Minkowski()
    """
    In this test I calculate the transition probability as a function of Ω, of an inertial observer
    in flat spacetime with a gaussian switching function.
    """
    λ = 1 
    σ = 1
    d = 5*σ   
    Ω = 1.0 
    ε_contour = 1e-2

    initial_τs, final_τs =  [-d, -d], [d, d]
    integrate(f::Function) = hcubature(f, initial_τs, final_τs, maxevals=500000 , rtol=rtol)[1]

    Ls = LinRange(σ/10, 2.5σ, 20)
    Ms = []
    for (i, L) in enumerate(Ls)
        XA, XB = initialize_trajs("flat", 0, L, b)
        W = deform(_Ws["flat"], ε_contour)
        Ws, D = initialize_distributions(W, XA, XB, false, ε_numeric_derivative)    
        χ(τ) = switching_funcs["gauss"](τ/σ)

        println("\rDoing L number $i: L = $L")
        m = get_m(D, λ, Ω, χ)
        M = integrate(m)
        push!(Ms, M)
    end

    M(L) = im*(λ^2)*σ/(4*√π*L)*exp(-(σ*Ω)^2 - L^2/(4*σ^2))*(erf(im*L/(2σ)) - 1)
    p = Plots.plot(Vec(Ls), [(abs ∘ M).(Ls), abs.(Ms)], labels=["theoretical" "numerical"], ylims=[-1e-6, max(abs.(Ms)...)*3/2])
    display(p)
    savefig(p, "plots\\test_plotsM_vs_L_Minkowski\\ε_contour_$(ε_contour)_σ=$σ.png")
end

function get_flat_concurrence()
    λ = 1.0
    σ = 1
    d = 5*σ    
    initial_τs, final_τs = [switching_function_center_A - d, switching_function_center_B - d], 
                           [switching_function_center_A + d, switching_function_center_B + d]
    ε_contour = 5e-3

    χ(τ) = switching_funcs["gauss"]((τ-switching_function_center_A)/σ)

    M_func(Ω, L) = im*(λ^2)*σ/(4*√π*L)*exp(-(σ*Ω)^2 - L^2/(4*σ^2))*(erf(im*L/(2σ)) - 1)
    P(Ω)         = λ^2/4π*(exp(-σ^2*Ω^2) - √π*σ*Ω*erfc(σ*Ω))
    C_func(Ω, L) = λ^2/4√π*σ/L*exp(-L^2/4σ^2) * (imag(exp(im*L*Ω) * erf(im*L/2σ + σ*Ω)) - sin(Ω*L))

    integrate = MemoizedIntegrator(initial_τs, final_τs, 500000, 1e-4)

    ΔLs = LinRange(0.5σ, 2σ, 10)
    Ωs  = LinRange(-3/σ, 3/σ, 12)
    Cs = zeros(length(Ωs), length(ΔLs))
    Cs_th = zeros(length(Ωs), length(ΔLs))
    for (i, Ω) in tqdm(enumerate(Ωs))
        for (j, ΔL) in tqdm(enumerate(ΔLs))
            XA, XB = initialize_trajs("flat", 0, ΔL, b)
            W = deform(_Ws["flat"], ε_contour)
            Ws, D = initialize_distributions(W, XA, XB, false, ε_numeric_derivative)    

            m, ls = get_m(D, λ, Ω, χ), get_ls(Ws, λ, Ω, χ)
            M, Ls = integrate(m)     , map_dict(integrate, ls)

            ρ_th = [         1 - 2P(Ω)                  0                0    M_func(Ω, ΔL);
                                    0                 P(Ω)    C_func(Ω, ΔL)               0;
                                    0  conj(C_func(Ω, ΔL))             P(Ω)               0;
                    conj(M_func(Ω, ΔL))                 0                0                0]

            ρ = get_ρ(M, Ls)
                    
            Cs[i,j]    = concurrence(ρ)
            Cs_th[i,j] = concurrence(ρ_th)
        end
    end

    p = Plots.plot(Plots.contourf(ΔLs, Ωs, Cs, ylabel="Ω", xlabel="ΔL"), Plots.contourf(ΔLs, Ωs, Cs_th, ylabel="Ω", xlabel="ΔL"), size=(3200,1800), linewidth=0, xtickfontsize=18, ytickfontsize=18)
    display(p)
    savefig(p, "plots\\test_plots\\flat_concurrence_heatmap.png")
    Cs, Cs_th
end

function plot_inertial_l()
    """
    In this test I calculate the transition probability as a function of Ω, of an inertial observer
    in flat spacetime with a gaussian switching function.
    """
    X = InertialTrajectory(0.0, 0.0, 0.0)
    W = DistributionWithTrajectories(_Ws["flat"], X, X)

    χ(τ) = switching_funcs["gauss"](τ/σ)
    df = deform_funcs["triangle"]
    Ω = 1
    l = get_l(W, λ, Ω, χ)
    l = complexify_l_or_m(l, df, distance_funcs["flat"], 1e-3)
    # l = complexify_l_or_m(l, 1e-5)
    
    _l(t) = l([t,0])
    _l_th(t) = λ^2*χ(t)*χ(0) * _W_flat_spacetime([t,0], [0,0])*exp(-im*Ω*t)

    ts = LinRange(-0.5, 0.5, 40)
    display(plot(ts, [abs ∘ _l, abs ∘ _l_th], ylims=[-4,4]))
end

function plot_inertial_m()
    """
    In this test I calculate the transition probability as a function of Ω, of an inertial observer
    in flat spacetime with a gaussian switching function.
    """
    L = 1.0
    XA, XB = InertialTrajectory(0.0, 0.0, 0.0), InertialTrajectory(L, 0.0, 0.0)
    W = DistributionWithTrajectories(_Ws["flat"], XA, XB)
    χ(τ) = switching_funcs["gauss"](τ/σ)

    df = deform_funcs["triangle"]
    Ω = 1
    m = get_m(W, λ, Ω, χ)
    m = complexify_l_or_m(m, df, distance_funcs["flat"], 1e-1)
    # m = complexify_l_or_m(l, 1e-5)
    
    _m(t) = m([t,0])
    _m_th(t) = -λ^2*χ(t)*χ(0) * _Ds["flat"]([t, 0.0, 0.0, 0.0], [0.0, L, 0.0, 0.0])*exp(im*Ω*t)
    ts = LinRange(-5, 5, 200)
    display(plot(ts, [abs ∘ _m, abs ∘ _m_th]))#, ylims=[-4,4]))
end

function accelerated_detector_tests()
    integrate = MemoizedIntegrator(initial_τs, final_τs, maxevals, rtol)
    χs = initialize_switching_funcs(switching_func_name, σ, switching_function_center_A, switching_function_center_B)  
    
    ρs = []
    Cs  ,Ns   = zeros(nΩ, nχ), zeros(nΩ, nχ)
    PBs ,PAs  = zeros(nΩ, nχ), zeros(nΩ)
    for (i, Ω) in tqdm(enumerate(Ωs))
        ΔLs = ΔLss[Ω]
        for (j, ΔL) in tqdm(enumerate(ΔLs))
            XA, XB = AcceleratedTrajectory(χ0, 0, 0, 0), AcceleratedTrajectory(χ0, ΔL, 0, 0)
            Ws, D  = initialize_distributions(_Ws[space_time], _Ds[space_time], XA, XB)
            if with_derivative_coupling Ws, D = add_crossed_derivatives(Ws, D, ε_numeric_derivative) end
    
            m, ls = get_m(D, λ, Ω, χs, ε_contour), get_ls(Ws, λ, Ω, χs, ε_contour)
            M, Ls = integrate(m)                 , map_dict(integrate, ls)
    
            ρ = get_ρ(M, Ls)
            push!(ρs, ρ)
            Cs[i,j], Ns[i,j]  = concurrence(ρ), negativity(ρ)
            PAs[i] , PBs[i,j] = real(ρ[2,2])  , real(ρ[3,3])
        end
    end
end

# get_P_Minkowski();
# Ωs, Ms_num, Ms_th = get_M_vs_Ω_Minkowski();
# get_M_vs_L_Minkowski();
Cs, Cs_th = get_flat_concurrence();
# plot_C_vs_L()
# plot_inertial_l();
# plot_inertial_m();