
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

function get_M_Minkowski()
    """
    In this test I calculate the transition probability as a function of Ω, of an inertial observer
    in flat spacetime with a gaussian switching function.
    """
    d = 3.5*σ    
    df = deform_funcs["cos2"]
    ε_contour = 1e-1
    L = 1.0

    initial_τs, final_τs =  [-d, -d], [d, d]
    integrate(f::Function) = hcubature(f, initial_τs, final_τs, maxevals=100000 , rtol=int_tol)[1]

    XA, XB = InertialTrajectory(0.0, 0.0, 0.0), InertialTrajectory(L, 0.0, 0.0)
    D = DistributionWithTrajectories(_Ds["flat"], XA, XB)
    χ(τ) = χs["gauss"](τ/σ)

    Ωs = LinRange(0, 5, 6)
    Ms = []
    for (i, Ω) in enumerate(Ωs)
        println("\rDoing Ω number $i")
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