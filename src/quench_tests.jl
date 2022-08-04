using HCubature, Plots, SpecialFunctions #, QuantumInformation
include("params.jl")
include("tools/D&W.jl")
include("tools/SwitchingFuncs.jl")
include("tools/DeformFuncs.jl")
include("tools/LM_getters.jl")

function get_P_Minkowski()
    """
    In this test I calculate the transition probability as a function of Ω, of an inertial observer
    in flat spacetime with a gaussian switching function.
    """
    d = 4
    initial_τs, final_τs =  [-d, -d], [d, d]
    integrate(f::Function) = hcubature(f, initial_τs, final_τs, maxevals=50000 , rtol=int_tol)[1]

    X = InertialTrajectory(0.0, 0.0, 0.0)
    W = DistributionWithTrajectories(_Ws["flat"], X, X)
    χ(τ) = χs["gauss"](τ/σ)
    df = deform_funcs[deform_func_name]

    Ωs = LinRange(0, 10, 20)
    P_Ms = []
    for (i, Ω) in enumerate(Ωs)
        println("\rDoing Ω number $i")
        l = get_l(W, λ, Ω, χ, distance_funcs["lorentz"], df, ε_contour)
        P_M = integrate(l)
        push!(P_Ms, P_M)
    end
    x = max(imag.(P_Ms)...)
    println("$x")

    P(Ω) = λ^2/4π*(exp(-σ^2*Ω^2) - √π*σ*Ω*erfc(σ*Ω))
    display(plot(Ωs, [P.(Ωs), real.(P_Ms)]))
end

function plot_inertial_l()
    """
    In this test I calculate the transition probability as a function of Ω, of an inertial observer
    in flat spacetime with a gaussian switching function.
    """
    d = 8

    X = InertialTrajectory(0.0, 0.0, 0.0)
    W = DistributionWithTrajectories(_Ws["flat"], X, X)
    χ(τ) = χs["gauss"](τ/σ)
    df = deform_funcs[deform_func_name]

    Ω = 1
    l = get_l(W, λ, Ω, χ, distance_funcs["lorentz"], df, ε_contour)
    _l(t) = l([t,0])
    _l_th(t) = λ^2*χ(t)*χ(0) * _W_flat_spacetime([t,0], [0,0])

    ts = LinRange(-2,2,100)
    display(plot(ts, [abs∘_l, abs∘_l_th], ylims=[0,2]))
end

function plot_accelerated_l()
    """
    In this test the two functions mostly coincide. The small difference in the curves I suspect 
    can be attributed to O(ε^2) corrections that birrel neglets when doing his manipulations
    """
    X = AcceleratedTrajectory(1.0, 0.0, 0.0)
    _W = _Ws["flat"]
    W = DistributionWithTrajectories(_W, X, X)

    d = 5
    iε = im*1e-1
    t = LinRange(0,d,500)

    l = real ∘ get_l(W, 1, 1, χs["one"])
    f(t) = l([t,0.0])
    g(t) = real(-cos(t)/(sinh(t/2 - 0.50*iε))^2)/16/π^2
    display(plot(t, [f, g], label=["l_mine" "l_birrel"]))
end

function get_P_minkowski()
    ε_W = 0.4
    d = 20*ε_W

    initial_τs, final_τs =  [-d, -d], [d, d]
    integrate(f::Function) = hcubature(f, initial_τs, final_τs, maxevals=100000 , rtol=int_tol)

    X = AcceleratedTrajectory(1.0, 0.0, 0.0)
    _W = _Ws["flat"]
    W = DistributionWithTrajectories(_W, X, X)

    probs, errs = [], []
    Ωs = LinRange(0.1,5,20)
    for Ω in Ωs
        l = real ∘ get_l(W, 1, Ω, χs["one"])
        x = integrate(l)./(2*d)
        println(x)
        push!(probs,x[1]); push!(errs,x[2])
    end

    over_boltz(Ω) = 0.68*Ω/(2π*(exp(2π*Ω/3) - 1))    
    boltz(Ω)      =      Ω/(2π*(exp(2π*Ω  ) - 1))
    display(plot(Ωs, [probs, boltz.(Ωs), over_boltz.(Ωs)], label=["Numeric" "Boltzmann" "over fitted Boltzmann"]))
    probs
end

get_P_Minkowski()
# plot_inertial_l()
# probs =  plot_probability_vs_Ω();
