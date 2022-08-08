using HCubature, Plots, SpecialFunctions
include("params.jl")
include("tools/2_D&W.jl")
include("tools/3_LM_getters.jl")
include("tools/SwitchingFuncs.jl")
include("tools/DeformFuncs.jl")

function get_P_Minkowski()
    """
    In this test I calculate the transition probability as a function of Ω, of an inertial observer
    in flat spacetime with a gaussian switching function.
    """
    d = 4*σ    
    df = deform_funcs["cos2"]
    ε_contour = 1e-1

    initial_τs, final_τs =  [-d, -d], [d, d]
    integrate(f::Function) = hcubature(f, initial_τs, final_τs, maxevals=100000 , rtol=int_tol)[1]

    X = InertialTrajectory(0.0, 0.0, 0.0)
    W = DistributionWithTrajectories(_Ws["flat"], X, X)
    χ(τ) = χs["gauss"](τ/σ)

    Ωs = LinRange(0, 3, 10)
    P_Ms = []
    for (i, Ω) in enumerate(Ωs)
        println("\rDoing Ω number $i")
        l = get_l(W, λ, Ω, χ)
        # l = complexify_l_or_m(l, ε_contour)
        l = complexify_l_or_m(l, df, distance_funcs["lorentz"], ε_contour)
        P_M = integrate(l)
        push!(P_Ms, P_M)
    end
    x = max(imag.(P_Ms)...)
    println("$x")

    P(Ω) = λ^2/4π*(exp(-σ^2*Ω^2) - √π*σ*Ω*erfc(σ*Ω))
    display(plot(Ωs, [P.(Ωs), real.(P_Ms)]))
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
    χ(τ) = χs["gauss"](τ/σ)

    Ωs = LinRange(0, 2.5, 30)
    Ms = []
    for (i, Ω) in enumerate(Ωs)
        println("\rDoing Ω number $i: Ω = $Ω")
        m = get_m(D, λ, Ω, χ)
        # m = complexify_l_or_m(m, ε_contour)
        m = complexify_l_or_m(m, df, distance_funcs["lorentz"], ε_contour)
        M = integrate(m)
        push!(Ms, M)
    end

    M(Ω) = im*(λ^2)*σ/(4*√π*L)*exp(-(σ*Ω)^2 - L^2/(4*σ^2))*(erf(im*L/(2σ)) - 1)
    display(plot(Ωs, [(abs ∘ M).(Ωs), abs.(Ms)], labels=["theoretical" "numerical"]))
    return Ωs, Ms, M.(Ωs)
end

function get_M_vs_L_Minkowski()
    """
    In this test I calculate the transition probability as a function of Ω, of an inertial observer
    in flat spacetime with a gaussian switching function.
    """
    d = 5*σ    
    df = deform_funcs["cos2"]
    ε_contour = 5e-3

    initial_τs, final_τs =  [-d, -d], [d, d]
    integrate(f::Function) = hcubature(f, initial_τs, final_τs, maxevals=100000 , rtol=int_tol)[1]

    Ls = LinRange(0.001, 7, 30)
    Ms = []
    for (i, L) in enumerate(Ls)
        XA, XB = InertialTrajectory(0.0, 0.0, 0.0), InertialTrajectory(L, 0.0, 0.0)
        D = DistributionWithTrajectories(_Ds["flat"], XA, XB)
        χ(τ) = χs["gauss"](τ/σ)

        println("\rDoing L number $i: L = $L")
        m = get_m(D, λ, Ω, χ)
        m = complexify_l_or_m(m, ε_contour)
        # m = complexify_l_or_m(m, df, distance_funcs["lorentz"], ε_contour)
        M = integrate(m)
        push!(Ms, M)
    end

    M(L) = im*(λ^2)*σ/(4*√π*L)*exp(-(σ*Ω)^2 - L^2/(4*σ^2))*(erf(im*L/(2σ)) - 1)
    display(plot(Ls, [(abs ∘ M).(Ls), abs.(Ms)], labels=["theoretical" "numerical"], ylims=[-1e-6, 2e-4]))
    return Ls, Ms, M.(Ls)
end

function plot_inertial_l()
    """
    In this test I calculate the transition probability as a function of Ω, of an inertial observer
    in flat spacetime with a gaussian switching function.
    """
    X = InertialTrajectory(0.0, 0.0, 0.0)
    W = DistributionWithTrajectories(_Ws["flat"], X, X)

    χ(τ) = χs["gauss"](τ/σ)
    df = deform_funcs["triangle"]
    Ω = 1
    l = get_l(W, λ, Ω, χ)
    l = complexify_l_or_m(l, df, distance_funcs["lorentz"], 1e-3)
    # l = complexify_l_or_m(l, 1e-5)
    
    _l(t) = l([t,0])
    _l_th(t) = λ^2*χ(t)*χ(0) * _W_flat_spacetime([t,0], [0,0])*exp(-im*Ω*t)

    ts = LinRange(-0.5, 0.5, 40)
    display(plot(ts, [abs ∘ _l, abs ∘ _l_th], ylims=[-4,4]))
end

function plot_inertial_m()
#     """
#     In this test I calculate the transition probability as a function of Ω, of an inertial observer
#     in flat spacetime with a gaussian switching function.
#     """
#     L = 1.0
#     XA, XB = InertialTrajectory(0.0, 0.0, 0.0), InertialTrajectory(L, 0.0, 0.0)
#     W = DistributionWithTrajectories(_Ws["flat"], XA, XB)
#     χ(τ) = χs["gauss"](τ/σ)

#     df = deform_funcs["triangle"]
#     Ω = 1
#     m = get_m(W, λ, Ω, χ)
#     m = complexify_l_or_m(m, df, distance_funcs["lorentz"], 1e-1)
#     # m = complexify_l_or_m(l, 1e-5)
    
#     _m(t) = m([t,0])
#     _m_th(t) = -λ^2*χ(t)*χ(0) * _Ds["flat"]([t, 0.0, 0.0, 0.0], [0.0, L, 0.0, 0.0])*exp(im*Ω*t)
#     ts = LinRange(-5, 5, 60)
#     display(plot(ts, [abs ∘ _m, abs ∘ _m_th]))#, ylims=[-4,4]))
end

get_P_Minkowski();
# Ωs, Ms_num, Ms_th = get_M_vs_Ω_Minkowski();
# Ls, Ms_num, Ms_th = get_M_vs_L_Minkowski();
# plot_inertial_l();
# plot_inertial_m();