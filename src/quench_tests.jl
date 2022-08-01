using HCubature, Plots #, QuantumInformation
include("tools/misc.jl")
include("tools/params.jl")
include("tools/D&W.jl")
include("tools/Trajectories.jl")
include("tools/SwitchingFuncs.jl")

function plot_inertial_l_and_get_its_integral()
    """
    In this test I calculate the transition probability of an inertial observer
    in flat spacetime. It should be as close to zero as possible. Also, I plot the
    absolute value of the integrand of the transition function next to it's theoretical 
    value.
    """
    d = 5
    initial_τs, final_τs =  [-d, -d], [d, d]
    integrate(f::Function) = hcubature(f, initial_τs, final_τs, maxevals=100000 , rtol=int_tol)

    X = InertialTrajectory(0.0, 0.0, 0.0)
    W = DistributionWithTrajectories(_Ws["flat"], X, X, dist_ε)

    l = get_l(W, λ, Ω, χs["one"])
    x = integrate(l)./(2*d)
    println(x)

    g(t) = l([t,0.0])
    iε = im*dist_ε
    f(t) = λ^2*exp(-im*Ω*t)/(4π^2*(t - iε)^2)
    
    t = LinRange(0,5*dist_ε,100)
    display(plot(t, [abs ∘ g, abs ∘ f]))
end

function plot_accelerated_l()
    """
    In this test the two functions mostly coincide. The small difference in the curves I suspect 
    can be attributed to O(ε^2) corrections that birrel neglets when doing his manipulations
    """
    X = AcceleratedTrajectory(1.0, 0.0, 0.0)
    _W = _Ws["flat"]
    W = DistributionWithTrajectories(_W, X, X, dist_ε)

    d = 5
    iε = im*1e-1
    t = LinRange(0,d,500)

    l = real ∘ get_l(W, 1, 1, χs["one"])
    f(t) = l([t,0.0])
    g(t) = real(-cos(t)/(sinh(t/2 - 0.50*iε))^2)/16/π^2
    display(plot(t, [f, g], label=["l_mine" "l_birrel"]))
end

function plot_probability_vs_Ω()
    dist_ε = 0.4
    d = 20*dist_ε
    initial_τs, final_τs =  [-d, -d], [d, d]
    integrate(f::Function) = hcubature(f, initial_τs, final_τs, maxevals=100000 , rtol=int_tol)

    X = AcceleratedTrajectory(1.0, 0.0, 0.0)
    _W = _Ws["flat"]
    W = DistributionWithTrajectories(_W, X, X, dist_ε)

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

# plot_inertial_l_and_get_its_integral()
# plot_accelerated_l()
probs =  plot_probability_vs_Ω();
