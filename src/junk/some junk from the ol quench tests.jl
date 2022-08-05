
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
