using HCubature, Plots, SpecialFunctions, ProgressBars
include("params.jl")
include("tools/3_LM_getters.jl")
include("tools/SwitchingFuncs.jl")
include("tools/DeformFuncs.jl")
include("tools/Memoized_Integrator.jl")

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
    integrate(f::Function) = hcubature(f, initial_τs, final_τs, maxevals=100000, rtol=int_tol)[1]

    X = InertialTrajectory(0.0, 0.0, 0.0)
    W = DistributionWithTrajectories(_Ws["flat"], X, X)
    χ(τ) = switching_funcs["gauss"](τ/σ)

    Ωs = LinRange(0, 2.5, 30)
    P_Ms = []
    for (i, Ω) in enumerate(Ωs)
        println("\rDoing Ω number $i")
        l = get_l(W, λ, Ω, χ, ε_contour)
        P_M = integrate(l)
        push!(P_Ms, P_M)
    end
    x = max(imag.(P_Ms)...)
    println("$x")

    P(Ω) = λ^2/4π*(exp(-σ^2*Ω^2) - √π*σ*Ω*erfc(σ*Ω))

    p = plot(Ωs, [P.(Ωs), real.(P_Ms), imag.(P_Ms)], labels=["Theoretical" "Real part of numerical" "Imaginary part of numerical"])
    display(p)
    savefig(p, raw"plots\P_vs_Ω_Minkowski.png")
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
    integrate(f) = hcubature(f, initial_τs, final_τs, maxevals=100000 , rtol=int_tol)[1]

    XA, XB = InertialTrajectory(0.0, 0.0, 0.0), InertialTrajectory(L, 0.0, 0.0)
    D = DistributionWithTrajectories(_Ds["flat"], XA, XB)
    χ(τ) = switching_funcs["gauss"](τ/σ)

    Ωs = LinRange(0, 2.5, 30)
    Ms = []
    for (i, Ω) in enumerate(Ωs)
        println("\rDoing Ω number $i: Ω = $Ω")
        m = get_m(D, λ, Ω, χ, ε_contour)
        M = integrate(m)
        push!(Ms, M)
    end

    M(Ω) = im*(λ^2)*σ/(4*√π*L)*exp(-(σ*Ω)^2 - L^2/(4*σ^2))*(erf(im*L/(2σ)) - 1)

    p = plot(Ωs, [(abs ∘ M).(Ωs), abs.(Ms)], labels=["theoretical" "numerical"])
    display(p)
    savefig(p, raw"plots\M_vs_Ω_Minkowski.png")
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
    integrate(f::Function) = hcubature(f, initial_τs, final_τs, maxevals=500000 , rtol=int_tol)[1]

    Ls = LinRange(σ/10, 2.5σ, 20)
    Ms = []
    for (i, L) in enumerate(Ls)
        XA, XB = InertialTrajectory(0.0, 0.0, 0.0), InertialTrajectory(L, 0.0, 0.0)
        D = DistributionWithTrajectories(_Ds["flat"], XA, XB)
        χ(τ) = switching_funcs["gauss"](τ/σ)

        println("\rDoing L number $i: L = $L")
        m = get_m(D, λ, Ω, χ, ε_contour)
        M = integrate(m)
        push!(Ms, M)
    end

    M(L) = im*(λ^2)*σ/(4*√π*L)*exp(-(σ*Ω)^2 - L^2/(4*σ^2))*(erf(im*L/(2σ)) - 1)
    p = plot(Ls, [(abs ∘ M).(Ls), abs.(Ms)], labels=["theoretical" "numerical"], ylims=[-1e-6, max(abs.(Ms)...)*3/2])
    display(p)
    savefig(p, "plots\\M_vs_L_Minkowski\\ε_contour_$(ε_contour)_σ=$σ.png")
end

function plot_C_vs_L()
    λ = 1.0 
    σ = 1.0
    d = 5*σ    
    Ω = 1
    ε_contour = 1e-2
    χ(τ) = switching_funcs["gauss"](τ/σ)

    initial_τs, final_τs =  [-d, -d], [d, d]
    integrate(f::Function) = hcubature(f, initial_τs, final_τs, maxevals=100000 , rtol=int_tol)[1]

    separations = LinRange(σ/10, 7σ, 20)
    LABs = []
    for separation in separations
        XA, XB = InertialTrajectory(0.0, 0.0, 0.0), InertialTrajectory(separation, 0.0, 0.0)
        Ws = initialize_Ws(_Ws["flat"], XA, XB)
        ls = get_ls(Ws, λ, Ω, χ, ε_contour)
        LAB = integrate(ls["AB"])
        push!(LABs, LAB)
    end

    C_func(Ω, L) = λ^2/4√π*σ/L*exp(-L^2/4σ^2)* (imag(exp(im*L*Ω) * erf(im*L/2σ + σ*Ω)) - sin(Ω*L))
    p = plot(separations, [abs.(LABs), abs.(C_func.(Ω, separations))])
    display(p)
    savefig(p, "plots\\C_vs_L_Minkowski\\ε_contour_$(ε_contour)_σ=$σ.png")
end

function get_concurrence()
    λ = 1.0
    σ = 1
    d = 0.5*σ    
    initial_τs, final_τs = [switching_function_center_A - d, switching_function_center_B - d], 
                           [switching_function_center_A + d, switching_function_center_B + d]
    ε_contour = 5e-3

    χ(τ) = switching_funcs["cos4"]((τ-switching_function_center_A)/σ)

    M_func(Ω, L) = im*(λ^2)*σ/(4*√π*L)*exp(-(σ*Ω)^2 - L^2/(4*σ^2))*(erf(im*L/(2σ)) - 1)
    P(Ω)         = λ^2/4π*(exp(-σ^2*Ω^2) - √π*σ*Ω*erfc(σ*Ω))
    C_func(Ω, L) = λ^2/4√π*σ/L*exp(-L^2/4σ^2) * (imag(exp(im*L*Ω) * erf(im*L/2σ + σ*Ω)) - sin(Ω*L))

    integrate = MemoizedIntegrator(initial_τs, final_τs, 500000, 1e-4)

    ΔLs = LinRange(0.5σ, 2σ, 5)
    Ωs  = LinRange(-3/σ, 3/σ, 5)
    Cs = zeros(length(Ωs), length(ΔLs))
    Cs_th = zeros(length(Ωs), length(ΔLs))
    for (i, Ω) in tqdm(enumerate(Ωs))
        for (j, ΔL) in tqdm(enumerate(ΔLs))
            XA, XB = initialize_trajs(space_time, 0.5, ΔL, 0.0)
            Ws, D = initialize_distributions(_Ws["flat"], _Ds["flat"], XA, XB)

            m, ls = get_m(D, λ, Ω, χ, ε_contour), get_ls(Ws, λ, Ω, χ, ε_contour)
            M, Ls = integrate(m)                , map_dict(integrate, ls)

            ρ_th = [         1 - 2P(Ω)                  0                0    M_func(Ω, ΔL);
                                    0                 P(Ω)    C_func(Ω, ΔL)               0;
                                    0  conj(C_func(Ω, ΔL))             P(Ω)               0;
                    conj(M_func(Ω, ΔL))                 0                0                0]

            ρ = get_ρ(M, Ls)
                    
            Cs[i,j]    = concurrence(ρ)
            Cs_th[i,j] = concurrence(ρ_th)
        end
    end

    p = plot(contourf(ΔLs, Ωs, Cs, ylabel="Ω", xlabel="ΔL"), contourf(ΔLs, Ωs, Cs_th, ylabel="Ω", xlabel="ΔL"), size=(3200,1800), linewidth=0, xtickfontsize=18, ytickfontsize=18)
    display(p)
    savefig(p, "plots\\old_plots\\flat_concurrence_heatmap.png")
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

# get_P_Minkowski();
# Ωs, Ms_num, Ms_th = get_M_vs_Ω_Minkowski();
# get_M_vs_L_Minkowski();
Cs, Cs_th = get_concurrence();
# plot_C_vs_L()
# plot_inertial_l();
# plot_inertial_m();